#!/usr/bin/env python
#
# This GIANT code is adapted from:
# Copyright (C) 2013 Radim Rehurek <me@radimrehurek.com>
# Licensed under the GNU LGPL v2.1 - http://www.gnu.org/licenses/lgpl.html

from __future__ import with_statement

import logging

from . import utils, matutils
from six.moves import xrange


logger = logging.getLogger('gensimmod.interfaces')


class CorpusABC(utils.SaveLoad):
    def __iter__(self):
        """
        Iterate over the corpus, yielding one document at a time.
        """
        raise NotImplementedError('cannot instantiate abstract base class')


    def save(self, *args, **kwargs):
        import warnings
        warnings.warn("corpus.save() stores only the (tiny) iteration object; "
            "to serialize the actual corpus content, use e.g. MmCorpus.serialize(corpus)")
        super(CorpusABC, self).save(*args, **kwargs)

    def __len__(self):
        """
        Return the number of documents in the corpus.

        This method is just the least common denominator and should really be
        overridden when possible.
        """
        raise NotImplementedError("must override __len__() before calling len(corpus)")
#        logger.warning("performing full corpus scan to determine its length; was this intended?")
#        return sum(1 for doc in self) # sum(empty generator) == 0, so this works even for an empty corpus

    @staticmethod
    def save_corpus(fname, corpus, id2word=None, metadata=False):
        """
        Save an existing `corpus` to disk.

        Some formats also support saving the dictionary (`feature_id->word` mapping),
        which can in this case be provided by the optional `id2word` parameter.

        >>> MmCorpus.save_corpus('file.mm', corpus)

        Some corpora also support an index of where each document begins, so
        that the documents on disk can be accessed in O(1) time (see the
        `corpora.IndexedCorpus` base class). In this case, `save_corpus` is automatically
        called internally by `serialize`, which does `save_corpus` plus saves the index
        at the same time, so you want to store the corpus with::

        >>> MmCorpus.serialize('file.mm', corpus) # stores index as well, allowing random access to individual documents

        Calling `serialize()` is preferred to calling `save_corpus()`.

        """
        raise NotImplementedError('cannot instantiate abstract base class')

        # example code:
        logger.info("converting corpus to ??? format: %s" % fname)
        with utils.smart_open(fname, 'wb') as fout:
            for doc in corpus: # iterate over the document stream
                fmt = str(doc) # format the document appropriately...
                fout.write(utils.to_utf8("%s\n" % fmt)) # serialize the formatted document to disk
#endclass CorpusABC


class TransformedCorpus(CorpusABC):
    def __init__(self, obj, corpus, chunksize=None):
        self.obj, self.corpus, self.chunksize = obj, corpus, chunksize
        self.metadata = False

    def __len__(self):
        return len(self.corpus)

    def __iter__(self):
        if self.chunksize:
            for chunk in utils.grouper(self.corpus, self.chunksize):
                for transformed in self.obj.__getitem__(chunk, chunksize=None):
                    yield transformed
        else:
            for doc in self.corpus:
                yield self.obj[doc]

    def __getitem__(self, docno):
        if hasattr(self.corpus, '__getitem__'):
           return self.obj[self.corpus[docno]]
        else:
            raise RuntimeError('Type {} does not support slicing.'.format(type(self.corpus)))
#endclass TransformedCorpus


class TransformationABC(utils.SaveLoad):
    def __getitem__(self, vec):
        """
        Transform vector from one vector space into another

        **or**

        Transform a whole corpus into another.
        """
        raise NotImplementedError('cannot instantiate abstract base class')


    def _apply(self, corpus, chunksize=None):
        """
        Apply the transformation to a whole corpus (as opposed to a single document)
        and return the result as another corpus.
        """
        return TransformedCorpus(self, corpus, chunksize)
#endclass TransformationABC


class SimilarityABC(utils.SaveLoad):
    def __init__(self, corpus):
        raise NotImplementedError("cannot instantiate Abstract Base Class")


    def get_similarities(self, doc):
        # (Sparse)MatrixSimilarity override this method so that they both use the
        # same  __getitem__ method, defined below
        raise NotImplementedError("cannot instantiate Abstract Base Class")


    def __getitem__(self, query):
        """Get similarities of document `query` to all documents in the corpus.

        **or**

        If `query` is a corpus (iterable of documents), return a matrix of similarities
        of all query documents vs. all corpus document. Using this type of batch
        query is more efficient than computing the similarities one document after
        another.
        """
        is_corpus, query = utils.is_corpus(query)
        if self.normalize:
            # self.normalize only works if the input is a plain gensim vector/corpus (as
            # advertised in the doc). in fact, input can be a numpy or scipy.sparse matrix
            # as well, but in that case assume tricks are happening and don't normalize
            # anything (self.normalize has no effect).
            if matutils.ismatrix(query):
                import warnings
                # warnings.warn("non-gensim input must already come normalized")
            else:
                if is_corpus:
                    query = [matutils.unitvec(v) for v in query]
                else:
                    query = matutils.unitvec(query)
        result = self.get_similarities(query)

        if self.num_best is None:
            return result

        # if the input query was a corpus (=more documents), compute the top-n
        # most similar for each document in turn
        if matutils.ismatrix(result):
            return [matutils.full2sparse_clipped(v, self.num_best) for v in result]
        else:
            # otherwise, return top-n of the single input document
            return matutils.full2sparse_clipped(result, self.num_best)


    def __iter__(self):
        """
        For each index document, compute cosine similarity against all other
        documents in the index and yield the result.
        """
        # turn off query normalization (vectors in the index are assumed to be already normalized)
        norm = self.normalize
        self.normalize = False

        # Try to compute similarities in bigger chunks of documents (not
        # one query = a single document after another). The point is, a
        # bigger query of N documents is faster than N small queries of one
        # document.
        #
        # After computing similarities of the bigger query in `self[chunk]`,
        # yield the resulting similarities one after another, so that it looks
        # exactly the same as if they had been computed with many small queries.
        try:
            chunking = self.chunksize > 1
        except AttributeError:
            # chunking not supported; fall back to the (slower) mode of 1 query=1 document
            chunking = False
        if chunking:
            # assumes `self.corpus` holds the index as a 2-d numpy array.
            # this is true for MatrixSimilarity and SparseMatrixSimilarity, but
            # may not be true for other (future) classes..?
            for chunk_start in xrange(0, self.index.shape[0], self.chunksize):
                # scipy.sparse doesn't allow slicing beyond real size of the matrix
                # (unlike numpy). so, clip the end of the chunk explicitly to make
                # scipy.sparse happy
                chunk_end = min(self.index.shape[0], chunk_start + self.chunksize)
                chunk = self.index[chunk_start : chunk_end]
                for sim in self[chunk]:
                    yield sim
        else:
            for doc in self.index:
                yield self[doc]

        # restore old normalization value
        self.normalize = norm
