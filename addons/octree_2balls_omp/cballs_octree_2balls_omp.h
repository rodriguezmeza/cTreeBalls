#ifndef _cballs_octree_2balls_omp_h
#define _cballs_octree_2balls_omp_h

    case OCTREE2BALLSMETHOD:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with octree two-ball 2PCF and LogMultipole 3PCF\n\n");
        if (cballs_opt_read_mask(cmd)) {
            bool cached_tree;
            int leaf_capacity;

            ifile = gd->iCatalogs[0];
            DO_BODY(p, bodytable[ifile],
                    bodytable[ifile] + gd->nbodyTable[ifile])
                Update(p) = TRUE;
            leaf_capacity = cballs_native_pair_leaf_capacity(
                cmd, gd, gd->nbodyTable[ifile]);
            cached_tree = !scanopt(cmd->options, "no-native-tree-cache")
                && octree_2balls_tree_cache_contains(
                    cmd, bodytable[ifile], gd->nbodyTable[ifile],
                    leaf_capacity);
            if (!cached_tree
                && MakeTree(cmd, gd, bodytable[ifile],
                            gd->nbodyTable[ifile], ifile) == FAILURE)
                return FAILURE;
            if (searchcalc_octree_2balls_omp(
                    cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                    ifile, ifile) == FAILURE)
                return FAILURE;
        } else {
            for (ifile = 0; ifile < gd->ninfiles; ifile++) {
                bool cached_tree;
                int leaf_capacity;

                DO_BODY(p, bodytable[ifile],
                        bodytable[ifile] + gd->nbodyTable[ifile])
                    Update(p) = TRUE;
                leaf_capacity = cballs_native_pair_leaf_capacity(
                    cmd, gd, gd->nbodyTable[ifile]);
                cached_tree = !scanopt(cmd->options, "no-native-tree-cache")
                    && octree_2balls_tree_cache_contains(
                        cmd, bodytable[ifile], gd->nbodyTable[ifile],
                        leaf_capacity);
                if (!cached_tree
                    && MakeTree(cmd, gd, bodytable[ifile],
                                gd->nbodyTable[ifile], ifile) == FAILURE)
                    return FAILURE;
            }
            if (searchcalc_octree_2balls_omp(
                    cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                    gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
                return FAILURE;
        }
        break;

#endif /* !_cballs_octree_2balls_omp_h */
