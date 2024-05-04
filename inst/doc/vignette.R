## ---- setup, echo=FALSE-------------------------------------------------------
GITHUB_README <- Sys.getenv("GITHUB_README") != ""
CAN_IGRAPH_PLOT <- requireNamespace("igraph", quietly=TRUE) && requireNamespace("ggplot2", quietly=TRUE)
knitr::opts_chunk$set(dpi=96,fig.width=6.5)

## ---- basic_usage, eval=FALSE-------------------------------------------------
#  results <- dist_search(x, y, max_distance = 2, nthreads = 1)

## ---- basic_plot, eval=!GITHUB_README && CAN_IGRAPH_PLOT, out.width=400-------
library(seqtrie)
tree <- RadixTree$new()
tree$insert(c("cargo", "cart", "carburetor", "carbuncle", "bar", "zebra"))
tree$erase("zebra")
# tree$graph requires igraph package
set.seed(1); tree$graph()

## ---- basic_plot_output, eval=GITHUB_README && CAN_IGRAPH_PLOT, echo=FALSE, message=FALSE, results='hide'----
#  library(seqtrie)
#  tree <- RadixTree$new()
#  tree$insert(c("cargo", "cart", "carburetor", "carbuncle", "bar", "zebra"))
#  tree$erase("zebra")
#  png("simple_tree.png", width = 400*1.5, height = 300*1.5, res = 96)
#  set.seed(1); tree$graph()
#  dev.off()

## ---- basic_plot_github, eval=GITHUB_README, echo=FALSE, results='asis'-------
#  cat('![](vignettes/simple_tree.png "simple_tree")')

## ---- small_cdr3_ex-----------------------------------------------------------
# 130,000 "CDR3" sequences
set.seed(1)
data(covid_cdr3)
covid_cdr3 <- sample(covid_cdr3, 1000)
tree <- RadixTree$new()
tree$insert(covid_cdr3)
# Full data: 1 min
results <- tree$search(covid_cdr3, max_distance=2, mode="levenshtein", nthreads=2)

# Alternatively, instead of using the RadixTree object directly, you can use the
# dist_search function, which is a wrapper around the RadixTree object.
results <- dist_search(covid_cdr3, covid_cdr3, max_distance=2)

# The output is a data.frame mapping query (search sequences)
# and target (sequences inserted into the tree).
dplyr::filter(results, query != target)

## ---- lv_search---------------------------------------------------------------
# Full data: several seconds
results <- tree$search(covid_cdr3, max_fraction=0.035, mode="levenshtein", nthreads=2)
# Full data: 1 minute
results <- tree$search(covid_cdr3, max_fraction=0.06, mode="levenshtein", nthreads=2)
# Full data: 15-20 minutes
results <- tree$search(covid_cdr3, max_fraction=0.15, mode="levenshtein", nthreads=2)

## ---- hm_search---------------------------------------------------------------
# Full data: 1 second
results <- tree$search(covid_cdr3, max_fraction=0.035, mode="hamming", nthreads=2)
# Full data: several seconds
results <- tree$search(covid_cdr3, max_fraction=0.06, mode="hamming", nthreads=2)
# Full data: 1.5 minutes
results <- tree$search(covid_cdr3, max_fraction=0.15, mode="hamming", nthreads=2)

## ---- anchored_search---------------------------------------------------------
tree <- RadixTree$new()
tree$insert("CARTON")
tree$insert("CAR")
tree$insert("CARBON")
tree$search("CART", max_distance = 0, mode = "anchored")

## ---- custom_search-----------------------------------------------------------
tree <- RadixTree$new()
tree$insert(covid_cdr3)
# define a custom distance matrix - generate_cost_matrix is a convienence function
# gap and gap_open can be defined directly in the cost_matrix or as search method parameters
cost_mat <- generate_cost_matrix("ACGT", match=0, mismatch=5, gap=2, gap_open=1)
print(cost_mat)
# Perform a search. "Mode" can be either global or anchored.
results <- tree$search(covid_cdr3, max_distance=8, cost_matrix=cost_mat, mode="global", nthreads=2)
dplyr::filter(results, query != target)

## ----radix_forest-------------------------------------------------------------
# RadixTree, full data: 45 seconds
tree <- RadixTree$new()
tree$insert(covid_cdr3)
results_tree <- tree$search(covid_cdr3, max_distance=2, mode="levenshtein", nthreads=2)
# RadixForest, full data: 19 seconds
frst <- RadixForest$new()
frst$insert(covid_cdr3)
results_frst <- frst$search(covid_cdr3, max_distance=2, mode="levenshtein", nthreads=2)
# The results are the same, but order is not guaranteed
identical(
  dplyr::arrange(results_tree, query, target),
  dplyr::arrange(results_frst, query, target) )

## ---- prefix_search-----------------------------------------------------------
tree <- RadixTree$new()
tree$insert(c("cargo", "cart", "carburetor", "carbuncle", "bar"))
tree$prefix_search("car")

