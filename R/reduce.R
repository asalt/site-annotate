library(purrr)
library(future)
library(furrr)

n_workers <- future::availableCores() - 1


reduce_sites <- function(datas, geneinfo, n_workers = NA, ...) {
  if (is.na(n_workers)) n_workers <- future::availableCores() - 1
  future::plan(future::multisession, workers = n_workers)
  # from https://stackoverflow.com/a/69008715

  datas_reduced <- datas %>%
    purrr::map(~ {
      common_cols <- intersect(colnames(geneinfo), colnames(.x)) # Check against the first dataframe or ensure it's generalizable
      .data_chunks <- .x %>%
        group_nest(fifteenmer, !!!rlang::syms(common_cols), .key = "grouped_data") %>%
        dplyr::mutate(.worker_id = sample(1:n_workers, replace = T, size = nrow(.))) %>%
        dplyr::group_split(.worker_id, .keep = F)
      .data_chunks %>%
        furrr::future_map_dfr(
          function(.data) {
            tidyr::unnest(.data, grouped_data) %>%
              group_by(fifteenmer, !!!rlang::syms(common_cols), ) %>%
              summarise(
                out = across(matches("TMT_.*intensity$"), .fns = \(x) sum(x, na.rm = TRUE)),
                .groups = "drop"
              ) %>%
              ungroup() %>%
              unnest(out)
          }
        )
    })
  return(datas_reduced)
}
