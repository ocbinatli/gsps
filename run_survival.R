library(dplyr)
library(survival)
library(survminer)

data_path <- "./data"
result_path <- "../results"
survgraph_path <- "../survgraph"

pathways <- c("hallmark", "pid", "kegg")

cohorts <- c("ACC",  "BLCA", "BRCA", "CESC", "CHOL",
              "COAD", "DLBC", "ESCA", "GBM",  "HNSC",
              "KICH", "KIRC", "KIRP", "LAML", "LGG",
              "LIHC", "LUAD", "LUSC", "MESO", "OV",
              "PAAD", "PCPG", "PRAD", "READ", "SARC",
              "SKCM", "STAD", "TGCT", "THCA", "THYM",
              "UCEC", "UCS",  "UVM")

methods <- c("GSPS", "SNF", "CC", "iCluster")

cluster_colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6")

for(cohort in cohorts) {
  load(sprintf("%s/TCGA-%s.RData", data_path, cohort))
  for(pathway in pathways) {
    for(method in methods) {
      for(cluster_count in 2:10) {
        if (file.exists(sprintf("%s/surv_figure_cohort_%s_pathway_%s_method_%s_%s.Rdata", survgraph_path, cohort, pathway, method, cluster_count)) == FALSE) {
          if (file.exists(sprintf("%s/state_cohort_%s_pathway_%s_%s_%s.RData", result_path, cohort, pathway, method, cluster_count)) == TRUE) {
            print(sprintf("running %s cohort with %d clusters with ...", cohort, cluster_count))
            load(sprintf("%s/state_cohort_%s_pathway_%s_%s_%s.RData", result_path, cohort, pathway, method, cluster_count))
            common_patients <- intersect(rownames(TCGA$clinical)[which(is.na(TCGA$clinical$vital_status) == FALSE)], rownames(TCGA$mrna))
            Y <- TCGA$clinical[common_patients, c("vital_status", "days_to_death", "days_to_last_followup")]
            
            delta <- 1 * (Y$vital_status == "Alive")
            y <- matrix(NA, length(delta), 1)
            y[delta == 0] <- Y$days_to_death[delta == 0]
            y[delta == 1] <- Y$days_to_last_followup[delta == 1]
            
            valid_patients <- which(y > 0 & is.na(y) == FALSE)
            y <- y[valid_patients,,drop = FALSE]
            delta <- delta[valid_patients]
        
            survival_matrix <- data.frame(matrix(NA, length(delta), 3, dimnames = list(1:length(delta), c("status", "time", "cluster"))))
            survival_matrix[, "status"] <- 1 - delta
            survival_matrix[, "time"] <- y
          
            survival_matrix[, "cluster"] <- state$clustering[valid_patients]
            cl <- cluster_count

            surv.diff <- survdiff(Surv(time, status) ~ cluster, data = survival_matrix)
            p.val <- 1 - pchisq(surv.diff$chisq, length(surv.diff$n) - 1)
            surv.fit <- survfit(Surv(time, status) ~ cluster, data = survival_matrix)
            
            if (cl == 3) {
              g <- ggsurvplot(surv.fit, size = 1.8, font.main = 17, font.submain = 16, font.caption = 16, font.x = 15, font.y = 15, font.tickslab = 15, font.legend = 15, risk.table.fontsize = 6, data = survival_matrix, pval = TRUE, risk.table = "absolute", xlab = "Days", ylab = "Survival ratio", legend.title = "", tables.height = 0.5 * log10(cluster_count), palette = cluster_colors[1:cluster_count])
            } else {
              g <- ggsurvplot(surv.fit, size = 1.8, font.main = 17, font.submain = 16, font.caption = 16, font.x = 15, font.y = 15, font.tickslab = 15, font.legend = 15, risk.table.fontsize = 6, data = survival_matrix, pval = TRUE, risk.table = "absolute", xlab = "Days", ylab = "Survival ratio", legend.title = "", tables.height = 0.4 * log10(cluster_count), palette = cluster_colors[1:cluster_count])
            }
            save("g", file = sprintf("%s/surv_figure_cohort_%s_pathway_%s_method_%s_%s.RData", survgraph_path, cohort, pathway, method, cluster_count))
          }
        }
      }
    }
  }
}
