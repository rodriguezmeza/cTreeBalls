
#ifndef _tpc_io_gadget_02_h
#define _tpc_io_gadget_02_h

local int inputdata_gadget(struct cmdline_data* cmd, struct  global_data* gd)
{
//    global float * density;
    float * density;
    particle_data_ptr Ppointer[1];
//    int files;
    bodyptr p;
    int ip, k;
    real weight=1.0;
//    files=1;                               // Number of files per snapshot. Change it if needed

//    load_snapshot(cmd->infile,  files, Ppointer);
//    P=Ppointer[0];

    read_catalog_tpcf(cmd->infile, &NumPart, Ppointer);
    PP=Ppointer[0];

    GridScale=3; // How many times finer (in 1-dim) is the PM grid compared to the particle grid.
    NROW=(int) (0.5+powf(((float )NumPart),1./3.));
    NGRID=GridScale*NROW;
    Scale=2.*M_PI/BoxSize;
    printf("input_gadget: NROW,NGRID,Scale: %i %i %g\n",NROW,NGRID,Scale);

    density=malloc(NGRID*NGRID*NGRID*sizeof(float)); // PM density grid.
    PtoMesh(density, PP);

    DO_COORD(k)
    gd->Box[k] = BoxSize;

    cmd->nbody = NumPart;
    verb_print_debug(1, "\nAqui voy (0) :: cmd->nbody = %ld\n",cmd->nbody);
    fprintf(stdout,"Aqui voy (1):: Pos0 %g\n",PP[1].Pos[0]);
    fflush(stdout);

    int ifile=0;
    gd->nbodyTable[ifile] = cmd->nbody;

//    bodytab = (bodyptr) allocate(cmd->nbody * sizeof(body));
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

//    DO_BODY(p, bodytab, bodytab+cmd->nbody) {
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
//        ip = p-bodytab;
        ip = p-bodytable[ifile];
        DO_COORD(k)
            Pos(p)[k] = PP[ip].Pos[k];
        Kappa(p) = 1.0 + rcos(2.0*PI*Pos(p)[0]/(gd->Box[0]/20.0))
                        * rsin(2.0*PI*Pos(p)[1]/(gd->Box[1]/20.0));
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
        Type(p) = BODY;
        Weight(p) = weight;
//        Id(p) = p-bodytab+1;
        Id(p) = p-bodytable[ifile]+1;
    }

    verb_print_debug(1, "\nAqui voy (1) :: cmd->nbody = %ld\n",cmd->nbody);

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
//        Type(p) = BODY;
//        Weight(p) = weight;
//        Id(p) = p-bodytab+1;
        kavg += Kappa(p);
    }
    verb_print(cmd->verbose,
               "inputdata_gadget: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );

    return _SUCCESS_;
}

#endif	// ! _tpc_io_gadget_02_h
