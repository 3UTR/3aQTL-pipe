# perform fine mapping on one gene
args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
Lvalue <- as.integer(args[2])
sp_var <- as.numeric(args[3])
mPIP <- as.numeric(args[4])

setwd(dir)

cat('Running environment:',getwd(),'\nOptions:\n','Lvalue:',Lvalue,'\nsp_var:',sp_var,'\nmPIP:',mPIP,'\n')
genotype = '3aQTL.vcf'
phenotype = 'expr.phen'

# output file prefix
prefix = tools::file_path_sans_ext(genotype)

X = t(read.table(genotype, head=T, row.names=1, quote="'"))
# fill missing values in X with mean
# because susieR does not deal with missing data explicitly for now
for(i in 1:ncol(X)){
  X[is.na(X[,i]), i] <- mean(X[,i], na.rm = TRUE)
}
y = read.table(phenotype, head=F)[,-1]
# Adjust row names for phenotype data convention
#rownames(y) = gsub("-", ".", y[,1])
rownames(y) = y[,1]
# Obtain intersect of X and y data, and reorder X to match y ordering
x_idx = match(rownames(y), rownames(X))
y_idx = which(!is.na(x_idx))
x_idx = x_idx[!is.na(x_idx)]
X = X[x_idx, ]
y = y[y_idx,]
if (!all(rownames(X) == rownames(y))) stop("X and y rownames mismatch")
# Run SuSiE
res = susieR::susie(X, y[,2], L=Lvalue, scaled_prior_variance=sp_var)
# Visualize result
pdf(paste0(prefix, '.SuSiE.pdf'), width=10,height=5)
susieR::susie_plot(res, y = 'PIP')
dev.off()
# Format results focusing only on signals
res$var_names = colnames(X)
get_susie_output = function(unit, res, pip_cutoff = mPIP) {
        cs_id = cs_size = cs_purity = rep(NA, length(res$var_names))
        num_cs = length(res$sets$cs)
        for(id in 1:num_cs){
            idx = res$sets$cs[[id]]
            cs_id[idx] = names(res$sets$cs)[id]
            cs_size[idx] = length(res$sets$cs[[id]])
            cs_purity[idx] = res$sets$purity[id,1]
        }
        out = cbind.data.frame(rep(unit, length(res$var_names)),
                                res$var_names,
                                res$pip, cs_id, cs_size, cs_purity)
        colnames(out) = c("locus_id", "variant_id", "pip", "cs", "cs_size", "cs_purity")
        out[which(out[,3] >= pip_cutoff | !is.na(out[,4])), ]
    }
text_output = get_susie_output(genotype, res)
# Output to files
write.table(text_output, paste0(prefix, '.SuSiE.txt'), quote=FALSE, row.names=FALSE)
saveRDS(res, paste0(prefix, '.SuSiE.rds'))
