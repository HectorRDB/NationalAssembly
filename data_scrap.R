library(xml2)
library(downloader)
library(tidyverse)

# Get the list of all delegate ----

path <- "http://data.assemblee-nationale.fr/static/openData/repository/14/amo/deputes_senateurs_ministres_legislatures_XIV/AMO20_dep_sen_min_tous_mandats_et_organes_XIV.csv.zip"
download(path, destfile = "data/deputes.zip", mode = "wb")
unzip("data/deputes.zip", exdir = "data/")

delegate <- read_delim("data/acteurs.csv", delim = ";")
res <- NULL

for (i in 1:248){
  mandate <- str_c("mandats[1]/mandat[", i, "]/typeOrgane[1]")
  # Select people who's ith mandate is to be a depute
  df <- delegate %>% filter(UQ(as.name(mandate)) == "ASSEMBLEE")
  if(nrow(df) > 0){
    circo <- str_c("mandats[1]/mandat[", i, "]/election[1]/lieu[1]/numCirco[1]")
    dept <- str_c("mandats[1]/mandat[", i,
                  "]/election[1]/lieu[1]/numDepartement[1]")
    actor <- str_c("mandats[1]/mandat[", i, "]/acteurRef[1]")
    mandate <- str_c("mandats[1]/mandat[", i, "]/uid[1]")
    name <- "etatCivil[1]/ident[1]/prenom[1]"
    surname <- "etatCivil[1]/ident[1]/nom[1]"
    df <- df %>% select(UQ(as.name(circo)), UQ(as.name(dept)), UQ(as.name(actor)),
           UQ(as.name(name)), UQ(as.name(surname)), UQ(as.name(mandate)))
    # Build the data frame
    colnames(df) <- c("circo", "dept", "actor", "name", "surname", "mandate")
    df <- df %>% unite(identifiant, actor, mandate, sep = "") %>%
      mutate(circo = as.numeric(circo), dept = as.numeric(dept))
    # join
    res <- bind_rows(res, df)
  }
}
rm(circo, dept, actor, mandate, df, name, surname, i, path)
res <- res %>% mutate(iden = paste(dept, circo, sep = "-"), NbVote = 0,
                      NbYes = 0, NbNo = 0, NbAbst = 0, group = "")

# Get the list of votes and voting history ----

path="http://data.assemblee-nationale.fr/static/openData/repository/14/loi/scrutins/Scrutins_XIV.xml.zip"
download(path, destfile = "data/Scrutins_XIV.xml.zip", mode = "wb")
unzip("data/Scrutins_XIV.xml.zip", exdir = "data/")
liste <- xml_children(read_xml("data/Scrutins_XIV.xml"))
votes <- data.frame(Number = NA, Title = NA)

for (vote in liste){
  if(xml_text(xml_children(vote)[14]) == "DecompteNominatif"){
    number <- xml_text(xml_children(vote)[2])
    ventil <- xml_children(vote)[16]
    group <- xml_children(xml_children(xml_children(ventil)))
    title <- xml_text(xml_children(vote)[11])
    votes <- add_row(votes, Number = number, Title = title)
    temp <- data.frame(res$identifiant,NA)
    colnames(temp) <- c("identifiant", str_c("vote", number))
    for (b in group){
      intermediaire <- xml_children(xml_children(b))
      intermediaire <- xml_children(intermediaire[3])
      groupname <- xml_text(xml_children(b)[1])
      for (i in 1:4){
        deps <- data.frame(xml_text(
                              xml_children(intermediaire[i])))
        deps <- as.character(deps[,1])
        if(length(deps) != 0){
            if (i == 1){
              deps <- unlist(strsplit(deps, "MG"))
              deps <- unlist(strsplit(deps, "PSE"))
              deps <- unlist(strsplit(deps, "PAN"))
            }
            indices <- res$identifiant %in% deps
            indices_temp <- temp$identifiant %in% deps
            res[indices, "group"] <- paste(res[indices, "group"],
                                           as.character(groupname),
                                           sep = ",")

            if(i != 1){
              res[indices, "NbVote"] <- res[indices, ]$NbVote + 1
            }
            if (i == 1){
              temp[indices_temp, 2] <- "NotVoting"
            }
            if (i == 2){
              res[indices, "NbYes"] <- res[indices, "NbYes"] + 1
              temp[indices_temp, 2] <- "For"
            }
            if (i == 3){
              res[indices, "NbNo"] <- res[indices, "NbNo"] + 1
              temp[indices_temp, 2] <- "Against"
            }
            if (i == 4){
              res[indices, "NbAbst"] <- res[indices, "NbAbst"] + 1
              temp[indices_temp, 2] <- "Abstain"
            }
          }
      }
    }
    res <- data.frame(res, temp[, 2])
    colnames(res)[ncol(res)] <- paste0("vote", number)
  }
}
votes <- votes[-1, ]
rm(temp, ventil, groupname, intermediaire, number, group, b,
   deps, liste, i, title, vote, indices, indices_temp, path, delegate)

# Clean political groups to delete replicate ----
group <- res$group
group1 <- group2 <- group3 <- rep(0, length(group))
for(i in 1:length(group)){
  delegate_group <- unlist(str_split(group[i], ","))
  delegate_group <- delegate_group[!str_detect(delegate_group, '"')]
  delegate_group <- delegate_group[delegate_group != ""]
  delegate_group <- sort(unique(delegate_group), decreasing = T)
  group1[i] <- sort(delegate_group, decreasing = T)[1]
  group2[i] <- sort(delegate_group, decreasing = T)[2]
  group3[i] <- sort(delegate_group, decreasing = T)[3]
}
res$group1 <- group1; res$group2 <- group2; res$group3 <- group3
res <- res %>% filter(group1 != 0 & group1 != "" & (!is.na(group1))) %>%
               select(-group) %>%
               mutate(circo = as.integer(circo))

rm(i, group1, group2, group3, group, delegate_group)

# Clean group names ----
for(name in c("group1", "group2", "group3")){
  groupnames <- res[[name]]
  groupnames <- str_replace_all(groupnames, "PO656010", "UDI")
  groupnames <- str_replace_all(groupnames, "PO656018", "GDR")
  groupnames <- str_replace_all(groupnames, "PO656002", "SER")
  groupnames <- str_replace_all(groupnames, "PO707869", "LR")
  groupnames <- str_replace_all(groupnames, "PO656006", "LR")
  groupnames <- str_replace_all(groupnames, "PO684957", "LR")
  groupnames <- str_replace_all(groupnames, "PO713077", "SER")
  groupnames <- str_replace_all(groupnames, "PO656022", "RRPD")
  groupnames <- str_replace_all(groupnames, "PO645633", "Independent")
  groupnames <- str_replace_all(groupnames, "PO656014", "Greens")
  res[[name]] <- groupnames
}

# Clean the list of deputees appearing twice ----
doubles <- res %>% group_by(name, surname) %>% filter(n() > 1)
meta_double <- doubles %>% select(circo, dept, name, surname, identifiant,
                                  iden, group1, group2, group3) %>%
               group_by(name, surname) %>%
               summarise(circo = mean(circo, na.rm = T),
                         dept = mean(dept, na.rm = T),
                         iden = paste(circo, dept, sep = "-"),
                         identifiant = dplyr::first(identifiant),
                         group1 = unique(group1)[1],
                         group2 = unique(group2)[2],
                         group3 = unique(group2)[3])

cleaned <- meta_double %>% left_join(doubles %>%
                                     select(-one_of(colnames(meta_double)),
                                             -NbVote, -NbYes,
                                             -NbNo, -NbAbst)) %>%
                           gather(key = "votes", value = "voting", 9:653) %>%
                           filter(!is.na(voting)) %>%
                           distinct() %>%
                           spread(key = "votes", value = "voting") %>%
                           mutate(circo = as.integer(circo),
                                   dept = as.numeric(dept))


res <- res %>% anti_join(doubles) %>% full_join(cleaned)

write_csv(res, path = "data/voting_record.csv")
write_csv(votes, path = "data/votes.csv")
