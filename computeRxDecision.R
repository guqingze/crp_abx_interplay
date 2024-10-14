
# Function to compute prescription decisions based on antibiotic changes and routes
computeRxDecision <- function(data, abx_rank, route, num_abx) {
  
  # Shift variables to observe prior states of antibiotics, route, and number of antibiotics
  data[, prev_abx_rank := shift(abx_rank, type = "lag"), by = .(EpisodeID, ID)]
  data[, prev_route := shift(route, type = "lag"), by = .(EpisodeID, ID)]
  data[, prev_num_abx := shift(num_abx, type = "lag"), by = .(EpisodeID, ID)]

  # Calculate changes in antibiotic rank, route, and number of antibiotics
  data[, abx_rank_change := fcase(
    prev_abx_rank == abx_rank, 0,
    prev_abx_rank > abx_rank, -1,
    prev_abx_rank < abx_rank, 1
  )]
  
  data[, route_change := fcase(
    prev_route == route, 0,
    prev_route == "IV" & route == "Oral", -1,
    prev_route == "Oral" & route == "IV", 1
  )]
  
  data[, num_abx_change := fcase(
    prev_num_abx == num_abx, 0,
    prev_num_abx > num_abx, -1,
    prev_num_abx < num_abx, 1
  )]

  # Define the prescription decision based on changes
  data[, prescription_decision := fcase(
    abx_rank_change == 1, "Escalation",
    abx_rank_change == -1 & abx_rank != 0, "De-escalation",
    abx_rank_change == -1 & abx_rank == 0, "Stop",
    abx_rank_change == 0 & abx_rank == 1 & num_abx_change == -1, "De-escalation",
    abx_rank_change == 0 & abx_rank == 1 & num_abx_change == 1, "Escalation",
    abx_rank_change == 0 & abx_rank == 1 & num_abx_change == 0 & route_change == 1, "Escalation",
    abx_rank_change == 0 & abx_rank == 1 & num_abx_change == 0 & route_change == -1, "De-escalation",
    abx_rank_change == 0 & abx_rank == 1 & num_abx_change == 0 & route_change == 0, "Unchanged",
    default = NA_character_
  )]
  
  return(data)
}