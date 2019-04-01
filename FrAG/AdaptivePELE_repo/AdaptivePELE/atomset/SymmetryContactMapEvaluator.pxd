cdef class SymmetryContactMapEvaluator:
    cdef public list symmetries, proteinList, ligandList
    cdef public set symmetricAtoms
    cdef public dict symToRowMap
