# consolidate files

library(tidyverse)
library(here)
library(fs)

# load in annotated files
d_2343 = read_tsv(here('data/phase_changes/event_summaries_annotated/2343x1952.mvr.tsv')) %>% 
  mutate(cross = '2343x1952')
d_annotated = read_tsv(here('data/phase_changes/co_draws_all.tsv'))

# load in updated summaries - this is somewhat memory intensive
# ~1.5m rows in total
fnames = dir_ls(here('data/phase_changes/event_summaries'))

d_all = map(
  fnames,
  ~ read_tsv(., col_types = cols()) %>% 
    # filter(call == 'cross_over') %>% 
    mutate(
      indel_proximity = as.numeric(
        ifelse(indel_proximity == 'N/A', Inf, indel_proximity)
      )
    )
)

names(d_all) = str_extract(fnames, '[GB0-9]{4,5}x[0-9]{4}')
d_all_combined = bind_rows(d_all, .id = 'cross')

# semi joins to get all metrics for annotated crosses 
d_2343_updated = d_2343 %>% 
  select(check, comments, cross, read_name) %>% 
  left_join(d_all_combined, by = c('cross', 'read_name')) %>% 
  filter(check != 2) %>%  # only keep TPs and FPs, not borderline (8 out of 153, not incl. 10 missing)
  select(check, comments, cross, everything())

d_annotated_updated = d_annotated %>% 
  select(check, comments, cross, read_name) %>% 
  left_join(d_all_combined, by = c('cross', 'read_name')) %>% 
  filter(cross != '3071x2931', cross != '3071x3062')

# combine the two label sets
all_annotated = bind_rows(
  d_2343_updated,
  d_annotated_updated) %>% 
  filter(!is.na(chromosome)) # some bugged 2343x1952 reads that are now invalid to use

write_tsv(all_annotated, here('data/phase_changes/cos_all_annotated.tsv'))
