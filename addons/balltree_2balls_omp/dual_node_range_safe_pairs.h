/* Rare extreme-range fallback: O(N^2), O(radial bins) memory. Each MPI
 * rank has the replicated catalog; collective publication/failure is unchanged. */
static bool dual_node_requires_scaled_pairs(const fcfc_balltreeptr tree)
{
    for (INTEGER i=0;i<tree->npoint;i++) {
        const fcfc_ballpoint *p=&tree->packed_points[i];
        double values[2]={p->weight,p->kappa};
        for (int k=0;k<2;k++) if (values[k]!=0.0) {
            int exponent; frexp(values[k],&exponent);
            if (!isfinite(values[k]) || exponent>200 || exponent< -200) return TRUE;
        }
    }
    return FALSE;
}
static int dual_node_scaled_pairs(const dual_node_search_context *context,
                                  fcfc_balltreeptr a, fcfc_balltreeptr b, bool same)
{
    struct cmdline_data *cmd=context->cmd;
    struct global_data *gd=context->gd;
    size_t stride=(size_t)cmd->sizeHistN+1;
    cballs_scaled_sum *sums=calloc(2*stride,sizeof(*sums));
    int status=SUCCESS;
    if (sums==NULL || cballs_opt_smooth_pivot(cmd)) {
        snprintf(cmd->error_message,_ERRORMSGSIZE_,
                 "%s: numerical range: scaled pair accumulation requires scratch memory and unsmoothed fields",
                 DUAL_NODE_METHOD_NAME);
        free(sums); return FAILURE;
    }
    for (INTEGER i=0;i<a->npoint;i++) for (INTEGER j=same?i+1:0;j<b->npoint;j++) {
        const fcfc_ballpoint *p=&a->packed_points[i], *q=&b->packed_points[j];
        if (p->source && p->source==q->source) continue;
        const int bin=dual_node_bin_index_squared(context,
            dual_node_point_distance_squared(cmd,gd,p,q));
        if (bin<0) continue;
        double factors[4]={context->weighted?p->weight:1.0,
                           context->weighted?q->weight:1.0,p->kappa,q->kappa};
        cballs_scaled_add(&sums[bin],factors,2);
        cballs_scaled_add(&sums[stride+bin],factors,4);
        gd->histNN[bin]+=1.0; gd->histNNSubXi2pcf[bin]+=1.0; gd->nbbcalc++;
    }
    for (int bin=1;bin<=cmd->sizeHistN;bin++) {
        gd->histXi2pcf[bin]=cballs_scaled_ratio(sums[stride+bin],sums[bin]);
        if (!isfinite(gd->histXi2pcf[bin])) {
            snprintf(cmd->error_message,_ERRORMSGSIZE_,
                "%s: numerical range: normalized pair bin %d is not representable",DUAL_NODE_METHOD_NAME,bin);
            status=FAILURE; break;
        }
    }
    free(sums);
    if (status==SUCCESS) {
        gd->histogram_products |= CBALLS_PRODUCT_NN | CBALLS_PRODUCT_XI;
        if (cballs_opt_compute_histn(cmd) && cballs_opt_and_cf(cmd)) {
            for(int bin=1;bin<=cmd->sizeHistN;bin++) gd->histNN[bin]*=2.0;
            status=search_compute_HistN(cmd,gd,a->npoint);
        }
    }
    return status;
}
