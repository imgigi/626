# =============================================================================
# fit_attractiveness.R -- MCI's huff.fit() on the In-class 8 interaction matrix
#
# huff.fit() takes tens of minutes on 617 villages x 28 stores, because each
# iteration calls huff.attrac(), which recomputes the Huff shares for the whole
# matrix several times per store. It is therefore run once, here, and its
# result is saved to data/rds/huff_fit.rds; in-class-ex8.qmd reads that file
# and compares it with the direct calibration loop. The inputs are built
# exactly as in the document. Run from in-class-ex/in-class-ex8/.
# =============================================================================

pacman::p_load(sf, tidyverse, MCI)

stores <- st_read(dsn = "data/geospatial", layer = "STORES_TC", quiet = TRUE)
village <- st_read(dsn = "data/geospatial", layer = "VILLAGE_TC", quiet = TRUE)
pop2019 <- read_csv("data/aspatial/pop2019.csv", show_col_types = FALSE)
store_sales <- read_csv("data/aspatial/Store_sales_data.csv",
                        show_col_types = FALSE)

store_sales_aggre <- store_sales %>%
  group_by(`Store Code`) %>%
  summarise(Total_Bill = sum(Bills, na.rm = TRUE),
            Total_Sales = sum(Amount, na.rm = TRUE))

stores <- stores %>%
  left_join(store_sales_aggre, by = join_by(Store_CD == `Store Code`)) %>%
  mutate(Store_size = 100)

pop2019 <- pop2019 %>% mutate(V_ID = str_remove_all(V_ID, "-"))
village <- village %>%
  left_join(pop2019, by = join_by(VILLCODE == V_ID)) %>%
  filter(!is.na(P_CNT))

village_center <- st_point_on_surface(village)
distmat <- st_distance(village_center, stores)
rownames(distmat) <- village$VILLCODE
colnames(distmat) <- stores$Store_CD

ijMatrix <- distmat %>%
  as.data.frame() %>%
  rownames_to_column(var = "VILLCODE") %>%
  pivot_longer(cols = -VILLCODE, names_to = "Store_CD",
               values_to = "distance") %>%
  mutate(distance = as.numeric(distance)) %>%
  left_join(stores %>% st_drop_geometry() %>% select(Store_CD, Store_size),
            by = "Store_CD") %>%
  relocate(distance, .after = last_col())

village_market <- village %>%
  st_drop_geometry() %>%
  transmute(VILLCODE,
            demand = P_CNT * sum(stores$Total_Bill) / sum(P_CNT))

store_obs <- stores %>%
  st_drop_geometry() %>%
  select(Store_CD, Total_Bill)

# The same lambda the document estimates with huff.lambda().
dist_lambda <- huff.lambda(ijMatrix, "VILLCODE", "Store_CD", "Store_size", "distance",
                           gamma = 1, atype = "pow",
                           lambda_startv = -0.5, lambda_endv = -3, dtype = "pow",
                           localmarket_dataset = village_market,
                           origin_id = "VILLCODE", localmarket = "demand",
                           location_dataset = store_obs, location_id = "Store_CD",
                           location_total = "Total_Bill",
                           method = "bisection", iterations = 10)

t0 <- Sys.time()
huff_fit <- huff.fit(ijMatrix, "VILLCODE", "Store_CD", "Store_size", "distance",
                     lambda = dist_lambda$lambda,
                     localmarket_dataset = village_market,
                     origin_id = "VILLCODE", localmarket = "demand",
                     location_dataset = store_obs, location_id = "Store_CD",
                     location_total = "Total_Bill",
                     tolerance = 5, iterations = 3, output = "total")
runtime <- difftime(Sys.time(), t0, units = "mins")

write_rds(list(total = huff_fit, runtime_min = as.numeric(runtime)),
          "data/rds/huff_fit.rds")
cat("huff.fit finished in", round(as.numeric(runtime), 1), "minutes\n")
