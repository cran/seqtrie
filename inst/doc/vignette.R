## ----setup, echo=FALSE--------------------------------------------------------
GITHUB_README <- Sys.getenv("GITHUB_README") != ""
knitr::opts_chunk$set(dpi=96,fig.width=6.5)
library(seqtrie)

## ----basic_usage, eval=FALSE--------------------------------------------------
# data(covid_cdr3)
# results <- dist_search(covid_cdr3, max_distance = 3,
#                        nthreads = 8, tree_class = "StarTree")

## ----tree_benchmark, eval=!GITHUB_README, echo=FALSE, out.width="100%", fig.cap="Global edit-distance self-join benchmark with max_distance = 3 and nthreads = 8."----
knitr::include_graphics("vignette_benchmark.png")

## ----tree_benchmark_github, eval=GITHUB_README, echo=FALSE, results='asis'----
# cat('![](vignettes/vignette_benchmark.png "vignette_benchmark")')

## ----basic_plot, eval=FALSE---------------------------------------------------
# tree <- radix_tree()
# insert(tree, c("cargo", "cart", "carburetor", "carbuncle", "bar", "zebra"))
# erase(tree, "zebra")
# # plot_tree requires igraph and ggplot2
# set.seed(1); plot_tree(tree)

## ----basic_plot_static, eval=!GITHUB_README, echo=FALSE, out.width=400--------
knitr::include_graphics("simple_tree.png")

## ----basic_plot_github, eval=GITHUB_README, echo=FALSE, results='asis'--------
# cat('![](vignettes/simple_tree.png "simple_tree")')

## ----cdr3_setup, echo=FALSE---------------------------------------------------
# 130,000 "CDR3" sequences
set.seed(1)
data(covid_cdr3)
covid_cdr3 <- sample(covid_cdr3, 1000)
tree <- radix_tree()
insert(tree, covid_cdr3)

## ----hm_search----------------------------------------------------------------
results <- align_search(tree, covid_cdr3, max_fraction = 0.035,
                        mode = "hamming", nthreads = 2)
results <- align_search(tree, covid_cdr3, max_fraction = 0.06,
                        mode = "hamming", nthreads = 2)
results <- align_search(tree, covid_cdr3, max_fraction = 0.15,
                        mode = "hamming", nthreads = 2)

## ----anchored_search----------------------------------------------------------
tree <- radix_tree()
insert(tree, "CARTON")
insert(tree, "CAR")
insert(tree, "CARBON")
align_search(tree, "CART", max_distance = 0, mode = "anchored")

## ----custom_search------------------------------------------------------------
tree <- radix_tree()
insert(tree, covid_cdr3)

# Define a custom substitution matrix. Use generate_cost_matrix for convenience.
cost_mat <- generate_cost_matrix("ACGT", match = 0, mismatch = 5)
print(cost_mat)

# Set gap penalties via parameters (not in the matrix):
# - Linear gaps: set gap_cost only
# - Affine gaps: set both gap_cost and gap_open_cost

# Linear example
results_linear <- align_search(tree, covid_cdr3, max_distance = 8,
                               mode = "global",
                               cost_matrix = cost_mat,
                               gap_cost = 2,
                               nthreads = 2)

# Affine example
results_affine <- align_search(tree, covid_cdr3, max_distance = 8,
                               mode = "global",
                               cost_matrix = cost_mat,
                               gap_cost = 2,
                               gap_open_cost = 5,
                               nthreads = 2)

results_linear[results_linear$query != results_linear$target, , drop = FALSE]

## ----startree-----------------------------------------------------------------
st <- star_tree(c("ACGT", "ACGA", "AAAA", "AAAT"),
                max_distance = 1,
                mismatch_cost = 1,
                gap_cost = 1,
                nthreads = 2)
result(st)

# Search another query set using the same fixed costs and threshold.
align_search(st, c("ACGT", "AAAC"))

# The same path is available through dist_search().
dist_search(c("ACGT", "ACGA", "AAAA", "AAAT"),
            max_distance = 1,
            tree_class = "StarTree")

## ----anchored_startree--------------------------------------------------------
ast <- star_tree(c("ACGT", "ACG", "AAAA", "AA"),
                 max_distance = 1,
                 mode = "anchored",
                 mismatch_cost = 1,
                 gap_cost = 1,
                 nthreads = 2)
result(ast)

align_search(ast, c("ACGT", "AA"))

dist_search(c("ACGT", "ACG", "AAAA", "AA"),
            max_distance = 1,
            mode = "anchored",
            tree_class = "StarTree")

## ----hamming_startree---------------------------------------------------------
hst <- star_tree(c("ACGT", "ACGA", "TCGT", "ACG"),
                 max_distance = 1,
                 mode = "hamming",
                 nthreads = 2)
result(hst)

align_search(hst, c("ACGT", "TTGT"))

dist_search(c("ACGT", "ACGA", "TCGT", "ACG"),
            max_distance = 1,
            mode = "hamming",
            tree_class = "StarTree")

