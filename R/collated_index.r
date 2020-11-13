collated_index_fix <- function(data, s_sp, sindex_value = "SINDEX", bootID=NULL, boot_ind=NULL, glm_weights=TRUE, rm_zero=TRUE){

  data <- setDT(data)
  data = data[SPECIES == s_sp, ]
  y = as.numeric(as.character(unique(data$M_YEAR)))

  if(is.null(boot_ind)){
    site_id_ <- unique(data[, SITE_ID])
    boot_ind_ <- matrix(NA, nrow = 1, ncol = length(site_id_))
    boot_ind_[1, ] <- seq_len(length(site_id_))
    boot_ind <- list(site_id = site_id_, boot_ind = boot_ind_)
  } else {
    boot_ind = boot_ind
  }

  if(is.null(bootID)){
    bootID <- 0
  }else{
    booID <- bootID
  }

  if(bootID == 0){
    if(is.null(boot_ind)){
    a <- data.table(uID = as.character(paste0("s_",seq_len(data[, uniqueN(SITE_ID)]))), SITE_ID = as.character(unique(data[, SITE_ID])))
    }else{
    a <- data.table(uID = as.character(paste0("s_",seq_along(boot_ind$boot_ind[bootID+1, ]))), SITE_ID = as.character(boot_ind$site_id[seq_along(boot_ind$boot_ind[bootID+1, ])]))
    }
  } else {
    a <- data.table(uID = as.character(paste0("s_",seq_along(boot_ind$boot_ind[bootID, ]))), SITE_ID = as.character(boot_ind$site_id[boot_ind$boot_ind[bootID, ]]))
  }

  s_y <- data.table(expand.grid(a$uID, y[order(y)]))
  names(s_y) <- c("uID", "M_YEAR")
  s_y[ , uID:=as.character(uID)]
  s_y[ , M_YEAR:=as.numeric(as.character(M_YEAR))]
  data[ , SITE_ID:=as.character(SITE_ID)]
  data[ , M_YEAR:=as.numeric(as.character(M_YEAR))]

  setkey(a, uID); setkey(s_y, uID)
  s_y <- merge(s_y, a, all.x=TRUE)

  setkey(s_y, SITE_ID, M_YEAR); setkey(data, SITE_ID, M_YEAR)

  if(sindex_value != 'SINDEX'){
    data_boot <- merge(s_y, data[ , .(SITE_ID, M_YEAR, SINDEX, TOTAL_NM)][, (sindex_value) := data[, get(sindex_value)]])
    } else {
    data_boot <- merge(s_y, data[ , .(SITE_ID, M_YEAR, SINDEX, TOTAL_NM)])
    }

  sum_data <- data.table( BOOTi = bootID,
                          M_YEAR = y[order(y)],
                          NSITE = data_boot[, uniqueN(SITE_ID), by = M_YEAR][order(M_YEAR), V1],
                          NSITE_OBS = data_boot[, sum(SINDEX>0), by = M_YEAR][order(M_YEAR), V1])

  setkey(sum_data, M_YEAR)

  if(nrow(sum_data[NSITE_OBS > 0,]) != 0){

    if(isTRUE(rm_zero)){
      zero_site <- data_boot[ , V1:=sum(SINDEX>0), uID][V1==0, uID]
      data_boot <- data_boot[!uID %in% zero_site, ]
    }

    pred_data_boot <- expand.grid(unique(data_boot$uID), unique(data_boot$M_YEAR))
    names(pred_data_boot) <- c("uID", "M_YEAR")

    if (data_boot[, uniqueN(M_YEAR)] == 1 & data_boot[, uniqueN(SITE_ID)] == 1) warning(paste0(sindex_value, " available for only one site and one year: no collated index computed for ", s_sp, " - ", unique(data_boot$M_YEAR)), call. = FALSE)
    
    if (data_boot[, uniqueN(M_YEAR)] == 1) {
      mod_form <- as.formula(paste0("I(round(", sindex_value, ")) ~ factor(uID) - 1"))
    } else {
      mod_form <- as.formula(ifelse(data_boot[, uniqueN(SITE_ID)] > 1, paste0("I(round(", sindex_value, ")) ~ factor(uID) + factor(M_YEAR) - 1"), paste0("I(round(", sindex_value, ")) ~ factor(M_YEAR) - 1")))
    }

    if(!isTRUE(glm_weights)){
      data_boot[ , weights := 1]
    }else{
      ## make sure that column TOTAL_NM exists
      data_boot[ , weights := TOTAL_NM]
    }

    col_index <- try(speedglm::speedglm(mod_form, data = data_boot, family = poisson(), weights = data_boot$weights), silent = TRUE)
      if (class(col_index)[1] == "try-error") {
        col_index <- try(speedglm::speedglm(mod_form, data = data_boot, family = poisson(), weights = data_boot$weights, method = "qr"), silent = TRUE)
        if (class(col_index)[1] == "try-error") {
          col_index <- try(glm(mod_form, data = data_boot, family = poisson(), weights = data_boot$weights), silent = TRUE)
            if (class(col_index)[1] == "try-error") {
              res <- sum_data[, COL_INDEX := NA]
              return(list(col_index = res[ , .(BOOTi, M_YEAR, NSITE, NSITE_OBS, COL_INDEX)], site_id = a$SITE_ID ))
            }
        }
     }

    pred_data_boot$fitted <- predict(col_index, newdata = pred_data_boot, type = "response")

    co_index_res <- data.table(pred_data_boot)[order(M_YEAR), sum(fitted)/length(unique(data$SITE_ID)), by = M_YEAR]
    setnames(co_index_res, "V1", "COL_INDEX"); setkey(co_index_res, M_YEAR)
    res <- merge(sum_data, co_index_res, all.x = TRUE)[is.na(COL_INDEX), COL_INDEX := 0]

  } else {
    res <- sum_data[, COL_INDEX := 0]
  }

  return(list(col_index = res[ , .(BOOTi, M_YEAR, NSITE, NSITE_OBS, COL_INDEX)], site_id = a$SITE_ID ))
}