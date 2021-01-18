from typing import List
import hail as hl
import matplotlib.pyplot as plt

from .aggregation import Aggregation
from .plotting import plot_agg_hist
from .agg import agg_call_rate

# NEW CONVENTION: `mt.excluded_row`, `mt.excluded_col` are booleans indicating rows or
# columns excluded from analysis.

# NEW CONVENTION: Every "analysis component" accepts a `name` argument and it may only add one field
# with that name to the MatrixTable.

# NEW CONVENTION: An analysis component that wants to perform aggregations, should anntoate the
# currently excluded rows so it can appropriately filter its aggregation.

# NOTE: All-in methods like variant_qc and sample_qc are not allowed any more.

###############################################################################
## Analysis Components

def fvcr(aggregations: List[Aggregation],
         name: str,
         mt: hl.MatrixTable,
         min_call_rate: float):
    cr = agg_call_rate(mt)
    fails_filter = (cr < min_call_rate)
    newly_excluded = ~mt.exclude_row & fails_filter
    mt = mt.annotate_rows(
        exclude_row = mt.exclude_row | fails_filter,
        **{name: hl.struct(
            call_rate = cr,
            excluded_before = mt.exclude_row,
            newly_excluded = newly_excluded)})

    def aggregate(mt: hl.MatrixTable):
        return hl.struct(
            excluded_variants = hl.agg.filter(
                mt[name].newly_excluded,
                hl.agg.collect(mt.row_key)),
            call_rate_hist = hl.agg.filter(
                ~(mt[name].excluded_before),
                hl.agg.hist(mt[name].call_rate, 0.0, 1.0, 20)))

    def plot(data: hl.Struct, c, e):
        excluded_variants = data.excluded_variants
        n_excluded_variants = len(excluded_variants)
        plot_agg_hist(data.call_rate_hist)
        plt.title(name)
        plt.xlabel('call_rate')
        plt.ylabel('frequency')

    aggregations.append(Aggregation(plot, row_aggs=aggregate))
    return mt


def fscr(aggregations: List[Aggregation],
         name: str,
         mt: hl.MatrixTable,
         min_call_rate: float):
    cr = agg_call_rate(mt)
    fails_filter = (cr < min_call_rate)
    newly_excluded = ~mt.exclude_col & fails_filter
    mt = mt.annotate_cols(
        exclude_col = mt.exclude_col | fails_filter,
        **{name: hl.struct(
            call_rate = cr,
            excluded_before = mt.exclude_col,
            newly_excluded = newly_excluded)})

    def aggregate(mt: hl.MatrixTable):
        return hl.struct(
            excluded_samples = hl.agg.filter(
                mt[name].newly_excluded,
                hl.agg.collect(mt.col_key)),
            call_rate_hist = hl.agg.filter(
                ~(mt[name].excluded_before),
                hl.agg.hist(mt[name].call_rate, 0.0, 1.0, 20)))

    def plot(r, data: hl.Struct, e):
        excluded_samples = data.excluded_samples
        n_excluded_samples = len(excluded_samples)
        plot_agg_hist(data.call_rate_hist)
        plt.title(name)
        plt.xlabel('call_rate')
        plt.ylabel('frequency')

    aggregations.append(Aggregation(plot, col_aggs=aggregate))
    return mt


###############################################################################
## Set up test data

mt = hl.read_matrix_table(
    '/Users/dking/projects/hail-data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.mt')
mt = mt._filter_partitions(list(range(4)))
# add some missingness
mt = mt.annotate_rows(bad_row_p = hl.rand_unif(0, 1))
mt = mt.annotate_cols(bad_col_p = hl.rand_unif(0, 1))
mt = mt.annotate_entries(
    GT = hl.or_missing(hl.rand_bool(mt.bad_row_p) | hl.rand_bool(mt.bad_col_p),
                       mt.GT)
)

###############################################################################
## Set up exclude_row and exclude_col and aggs

mt = mt.annotate_rows(exclude_row = False)
mt = mt.annotate_cols(exclude_col = False)

aggs: List[Aggregation] = []

###############################################################################
## Pipeline:

mt = fvcr(aggs, 'variant call rate < .80', mt, 0.80)
mt = fscr(aggs, 'sample call rate < .80', mt, 0.80)
mt = fvcr(aggs, 'variant call rate < .90', mt, 0.90)
mt = fscr(aggs, 'sample call rate < .90', mt, 0.90)
mt = fvcr(aggs, 'variant call rate < .99', mt, 0.99)
mt = fscr(aggs, 'sample call rate < .99', mt, 0.99)
mt = fvcr(aggs, 'variant call rate < .999', mt, 0.999)
mt = fscr(aggs, 'sample call rate < .999', mt, 0.999)

mt.write('/tmp/foo.mt', overwrite=True)
mt = hl.read_matrix_table('/tmp/foo.mt')

###############################################################################
## Perform Aggregations

row_agg_results = mt.aggregate_rows(hl.tuple((x.row_aggs(mt) for x in aggs)))
col_agg_results = mt.aggregate_cols(hl.tuple((x.col_aggs(mt) for x in aggs)))
entry_agg_results = mt.aggregate_entries(hl.tuple((x.entry_aggs(mt) for x in aggs)))

##
###############################################################################

###############################################################################
## Run callbacks

callback_results = [
    agg.callback(row_result, col_result, entry_result)
    for (row_result, col_result, entry_result, agg)
    in zip(row_agg_results, col_agg_results, entry_agg_results, aggs)
]

##
###############################################################################

plt.show()
