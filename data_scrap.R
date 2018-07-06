# https://freakonometrics.hypotheses.org/50973
library(xml2)
library(downloader)
library(tidyverse)

# Get the list of all deputees ----

path <- "http://data.assemblee-nationale.fr/static/openData/repository/14/amo/deputes_senateurs_ministres_legislatures_XIV/AMO20_dep_sen_min_tous_mandats_et_organes_XIV.csv.zip"
download(path, destfile = "data/deputes.zip", mode = "wb")
unzip("data/deputes.zip", exdir = "data/")

deputes <- read_delim("data/acteurs.csv", delim = ";")
res <- NULL

for (i in 1:248){
  mandate <- str_c("mandats[1]/mandat[", i, "]/typeOrgane[1]")
  # Select people who's ith mandate is to be a depute
  df <- deputes %>% filter(UQ(as.name(mandate)) == "ASSEMBLEE")
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
                      NbYes = 0, NbNo = 0, NbAbst = 0, chamber = "")

# Get the list of votes and voting history ----

path="http://data.assemblee-nationale.fr/static/openData/repository/14/loi/scrutins/Scrutins_XIV.xml.zip"
download(path, destfile = "data/Scrutins_XIV.xml.zip", mode = "wb")
unzip("data/Scrutins_XIV.xml.zip", exdir = "data/")
liste <- xml_children(read_xml("data/Scrutins_XIV.xml"))
scrutins <- data.frame(Number = NA, Title = NA)

for (scrutin in liste){
  if(xml_text(xml_children(scrutin)[14]) == "DecompteNominatif"){
    number <- xml_text(xml_children(scrutin)[2])
    ventil <- xml_children(scrutin)[16]
    groupe <- xml_children(xml_children(xml_children(ventil)))
    title <- xml_text(xml_children(scrutin)[11])
    scrutins <- add_row(scrutins, Number = number, Title = title)
    temp <- data.frame(res$identifiant,NA)
    colnames(temp) <- c("identifiant", str_c("scrutin", number))
    for (b in groupe){
      intermediaire <- xml_children(xml_children(b))
      intermediaire <- xml_children(intermediaire[3])
      nomGroupe <- xml_text(xml_children(b)[1])
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
            res[indices, "chamber"] <- as.character(nomGroupe)
            
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
    colnames(res)[ncol(res)] <- paste("Scrutin", number)
  }
}
scrutins <- scrutins[-1, ]
res <- res %>% filter(chamber != 0)
rm(temp, ventil, nomGroupe, intermediaire, number, groupe, b,
   deps, liste, i, title, scrutin, indices, indices_temp, path)
write_csv(res, path = "data/voting_record.csv")
write_csv(scrutins, path = "data/scrutins.csv")