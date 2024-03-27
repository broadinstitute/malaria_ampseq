reverse_complement <- function(sequence) {
  complement <- c(A = "T", T = "A", C = "G", G = "C")
  paste(rev(complement[strsplit(toupper(sequence), "")[[1]]]), collapse = "")
}

search_forward = function(sequence, forward_barcode, forwardfound, forward_distance, five_prime_end, dist) {
  forward_barcode = as.character(forward_barcode)
  if (stringdist(substr(sequence, 1, nchar(forward_barcode)), forward_barcode, method = "hamming") <= dist) {
    forward_distance = stringdist(substr(sequence, 1, nchar(forward_barcode)), forward_barcode, method = "hamming")
    forwardfound = TRUE
    five_prime_end = substr(sequence, 1, nchar(forward_barcode))
  }

  return(list(forward_barcode = forward_barcode,
              forward_distance = forward_distance,
              five_prime_end = five_prime_end,
              forwardfound = forwardfound
  ))
}

search_reverse = function(sequence, reverse_barcode_original, reversefound, reverse_distance, three_prime_end, dist) {
  reverse_barcode = reverse_complement(as.character(reverse_barcode_original))
  if (stringdist(substr(sequence, nchar(sequence)-nchar(reverse_barcode)+1, nchar(sequence)), reverse_barcode, method = "hamming") <= dist && !reversefound) {
    reverse_distance = stringdist(substr(sequence, nchar(sequence)-nchar(reverse_barcode)+1, nchar(sequence)), reverse_barcode, method = "hamming")
    reversefound = TRUE
    three_prime_end = substr(sequence, nchar(sequence)-nchar(reverse_barcode)+1, nchar(sequence))
  }
  return(list(reverse_barcode = reverse_barcode_original,
              reverse_distance = reverse_distance,
              three_prime_end = three_prime_end,
              reversefound = reversefound
  ))
}

match_fun = function(sample_id, sequence, barcodes, dist) {
  forwardfound = FALSE
  reversefound = FALSE
  forward_distance = NA
  reverse_distance = NA
  five_prime_end = 'NULL'
  three_prime_end = 'NULL'
  forward_barcode = 'NULL'
  reverse_barcode = 'NULL'
  
  #Start by checking the expected barcodes
  forward = barcodes[which(barcodes$sample_id == sample_id), 'Forward']
  reverse = barcodes[which(barcodes$sample_id == sample_id), 'Reverse']
  match_tmp = search_forward(sequence, forward, forwardfound, forward_distance, five_prime_end, dist)
  if (match_tmp[[4]] == TRUE) {
    forward_barcode = match_tmp[[1]]
    forward_distance = match_tmp[[2]] 
    five_prime_end = match_tmp[[3]]
    forwardfound = match_tmp[[4]]
  }
  
  match_tmp = search_reverse(sequence, reverse, reversefound, reverse_distance, three_prime_end, dist)
  if (match_tmp[[4]] == TRUE) {
    reverse_barcode = match_tmp[[1]]
    reverse_distance = match_tmp[[2]] 
    three_prime_end = match_tmp[[3]]
    reversefound = match_tmp[[4]]
  }
  
  if (forwardfound && reversefound) {
    match = "Match"
    insert_size = nchar(sequence) - nchar(forward_barcode) - nchar(reverse_barcode)
    
  #If one or both of the expected barcodes are not found
    #First check it at the least one of the expected barcodes was fond
  } else if (forwardfound && !reversefound) {
    for (reverse in unique(barcodes$Reverse[which(barcodes$Reverse != reverse)])) {
      match_tmp = search_reverse(sequence, reverse, reversefound, reverse_distance, three_prime_end, dist)
      if (match_tmp[[4]] == TRUE) {
        reverse_barcode = match_tmp[[1]]
        reverse_distance = match_tmp[[2]] 
        three_prime_end = match_tmp[[3]]
        reversefound = match_tmp[[4]]
        match = "Forward_Match_Reverse_Mismatch"
        insert_size = nchar(sequence) - nchar(forward_barcode) - nchar(reverse_barcode)
        break
      } else {
        match = "Forward_Match_Reverse_Missing"
        insert_size = nchar(sequence) - nchar(forward_barcode) 
      }
    }  

  } else if (!forwardfound && reversefound) {
    for (forward in unique(barcodes$Forward[which(barcodes$Forward != forward)])) {
      match_tmp = search_forward(sequence, forward, forwardfound, forward_distance, five_prime_end, dist)
      if (match_tmp[[4]] == TRUE) {
        forward_barcode = match_tmp[[1]]
        forward_distance = match_tmp[[2]] 
        five_prime_end = match_tmp[[3]]
        forwardfound = match_tmp[[4]]
        match = "Forward_Mismatch_Reverse_Match"
        insert_size = nchar(sequence) - nchar(forward_barcode) - nchar(reverse_barcode)
        break
      } else{
        match = "Forward_Missing_Reverse_Match"
        insert_size = nchar(sequence) - nchar(reverse_barcode)
      }
    }  

    #If no barcodes were found, check all the other barcodes not assigned to the sample
  } else if (!forwardfound && !reversefound) {
    for (forward in unique(barcodes$Forward[which(barcodes$Forward != forward)])) {
      match_tmp = search_forward(sequence, forward, forwardfound, forward_distance, five_prime_end, dist)
      if (match_tmp[[4]] == TRUE) {
        forward_barcode = match_tmp[[1]]
        forward_distance = match_tmp[[2]] 
        five_prime_end = match_tmp[[3]]
        forwardfound = match_tmp[[4]]
        break
      }
    }  

    for (reverse in unique(barcodes$Reverse[which(barcodes$Reverse != reverse)])) {
      match_tmp = search_reverse(sequence, reverse, reversefound, reverse_distance, three_prime_end, dist)
      if (match_tmp[[4]] == TRUE) {
        reverse_barcode = match_tmp[[1]]
        reverse_distance = match_tmp[[2]] 
        three_prime_end = match_tmp[[3]]
        reversefound = match_tmp[[4]]
        break
      }
    }  

    if (forwardfound && reversefound) {
      match = "Mismatch" #Barcodes found, but they belong to a different well
      insert_size = nchar(sequence) - nchar(forward_barcode) - nchar(reverse_barcode)
    } else if (forwardfound && !reversefound) {
      match = "Forward_Mismatch_Reverse_Missing" #Forward found, but it belongs to a different well
      insert_size = nchar(sequence) - nchar(forward_barcode) 
      
    } else if (!forwardfound && reversefound) {
      match = "Reverse_Mismatch_Forward_Missing" #Reverse found, but they belong to a different well
      insert_size = nchar(sequence) - nchar(reverse_barcode)
      
    } else if (!forwardfound && !reversefound) {
      match = "No_barcodes" #No barcodes found 
      insert_size = nchar(sequence)
    }
    
  }
  
  return(list(forward_barcode = forward_barcode,
              reverse_barcode = reverse_barcode,
              forward_distance = forward_distance,
              reverse_distance = reverse_distance,
              five_prime_end = five_prime_end,
              three_prime_end = three_prime_end,
              match = match,
              insert_size = insert_size
  ))
}
