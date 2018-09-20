# Function to select region-specific sample IDs (independent of donors)
sample.ids <- function(roi, hemisphere = FALSE){
  sapply(roi, function(r){
    row <- match(r, ontology$name)
    id <- ontology$id[row]
    rows <- grep(id, ontology$structure_id_path)
    if (hemisphere) {
      hemisphere <- ontology$hemisphere[rows]
      rows_left <- rows[hemisphere == "L"]
      rows_right <- rows[hemisphere == "R"]
      list(left = ontology$id[rows_left], right = ontology$id[rows_right])
    }
    else ontology$id[rows]
  }, simplify = FALSE)
}