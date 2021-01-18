def agg_call_rate(mt: hl.MatrixTable):
    # DOES NOT HANDLE filter_entries CORRECTLY!
    n_called = hl.agg.count_where(hl.is_defined(mt['GT']))

    return hl.agg.filter(
        ~(mt.exclude_row | mt.exclude_col),
        n_called / hl.agg.count())
