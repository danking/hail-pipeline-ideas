class Aggregation:
    '''An Aggregation consists of four pieces:
- callback: a function which takes the aggregated rows, cols, and entries and three arguments
- row_aggs: None or a function taking a MatrixTable and returning an hl.struct of aggregations of
            row fields
- col_aggs: None or ... of column fields
- entry_aggs: None or ... of entry fields
'''

    def __init__(self,
                 callback: Callable[[hl.Struct, hl.Struct, hl.Struct], Any],
                 row_aggs: Optional[Callable[[hl.MatrixTable], Any]] = None,
                 col_aggs: Optional[Callable[[hl.MatrixTable], Any]] = None,
                 entry_aggs: Optional[Callable[[hl.MatrixTable], Any]] = None):
        if row_aggs is None:
            row_aggs = lambda mt: hl.struct()
        self.row_aggs = row_aggs

        if col_aggs is None:
            col_aggs = lambda mt: hl.struct()
        self.col_aggs = col_aggs

        if entry_aggs is None:
            entry_aggs = lambda mt: hl.struct()
        self.entry_aggs = entry_aggs

        self.callback = callback


