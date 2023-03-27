$# LS-DYNA Keyword file created by LS-PrePost(R) V4.7.17 - 08Jul2020
$# Created on Sep-15-2022 (18:47:29)
*KEYWORD
$ =============
$ CONTROL cards
$ =============
$
*CONTROL_ACCURACY
$$     OSU       INN    PIDOSU
$#     osu       inn    pidosu      iacc    
         0         4         0         0
*CONTROL_BULK_VISCOSITY
$$      Q2        Q1      TYPE
$#      q1        q2      type     btype    tstype      
       1.5      0.06        -1         0         0
*CONTROL_CONTACT
$#  slsfac    rwpnal    islchk    shlthk    penopt    thkchg     orien    enmass
       0.1       0.0         2         2         1         0         2         0
$#  usrstr    usrfrc     nsbcs    interm     xpene     ssthk      ecdt   tiedprj
         0         0        10         0       4.0         1         0         0
$#   sfric     dfric       edc       vfc        th     th_sf    pen_sf      
       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#  ignore    frceng   skiprwg    outseg   spotstp   spotdel   spothin       
         1         0         0         0         0         0       0.1
$#    isym    nserod    rwgaps    rwgdth     rwksf      icov    swradf    ithoff
         0         0         1       0.0       1.0         0       0.0         0
$#  shledg    pstiff    ithcnt    tdcnof     ftall    unused    shltrw    igactc
         0         0         0         0         0                 0.0         0
*CONTROL_ENERGY
$$    HGEN      RWEN    SLNTEN     RYLEN
$#    hgen      rwen    slnten     rylen     irgen     
         2         2         2         1         2
*CONTROL_HOURGLASS
$$     IHQ        QH
$#     ihq        qh  
         8      0.05
*CONTROL_OUTPUT
$$   NPOPT    NEECHO    NREFUP    IACCOP     OPIFS    IPNINT    IKEDIT    IFLUSH
$#   npopt    neecho    nrefup    iaccop     opifs    ipnint    ikedit    iflush
         1         3         0         1       0.0         0       100         0
$#   iprtf    ierode   tet10s8    msgmax    ipcurv      gmdt   ip1dblt      eocs
         0         0         2        50         0       0.0         0         0
$#   tolev    newleg    frfreq     minfo    solsig    msgflg    cdetol      
         2         0         1         0         0         0      10.0
$# phschng    demden   icrfile   spc2bnd    penout    shlsig   hisnout    engout
         0         0         0         0         0         0         0         0
$#    insf    isolsf      ibsf      issf    mlkbag      
         0         0         0         0         0
*CONTROL_SHELL
$$  WRPANG     ESORT     IRNXX    ISTUPD    THEORY       BWC     MITER      PROJ
$#  wrpang     esort     irnxx    istupd    theory       bwc     miter      proj
      20.0         1        -1         4         2         1         1         1
$$ ROTASCL    INTGRD    LAMSHT    CSTYP6    TSHELL    NFAIL1    NFAIL4
$# rotascl    intgrd    lamsht    cstyp6    thshel      
       1.0         0         0         1         0
$# psstupd   sidt4tu     cntco    itsflg    irquad    w-mode   stretch      icrq
         0         0         0         0         2       0.0       0.0         0
$#  nfail1    nfail4   psnfail    keepcs     delfr   drcpsid    drcprm   intperr
         1         1   6600000         0         0         0       0.0         0
*CONTROL_SOLID
$$   ESORT    FMATRX   NIPTETS    SWLOCL
$#   esort   fmatrix   niptets    swlocl    psfail   t10jtol    icohed    tet13k
         1         0         0         1         0       0.0         0         0
$#   pm1     pm2     pm3     pm4     pm5     pm6     pm7     pm8     pm9    pm10
       0       0       0       0       0       0       0       0       0       0
*CONTROL_SOLUTION
$#    soln       nlq     isnan     lcint     lcacc     ncdcf     
         0         0         1         0         0         1
$*CONTROL_TIMESTEP
$$  ENDTIM    ENDCYC     DTMIN    ENDENG    ENDMAS
$       0.0       0.9         0       0.0-5.0000E-4         0         0         0
$       0.0         0         0                           0.0                   0
$
$
$ ==============
$ DATABASE cards
$ ==============
$
*DATABASE_ABSTAT
$#      dt    binary      lcur     ioopt     
       0.1         2         0         1
*DATABASE_DEFORC
$#      dt    binary      lcur     ioopt     
       0.1         2         0         1
*DATABASE_ELOUT
$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
       0.1         2         0         1         0         0         0         0
*DATABASE_GLSTAT
$#      dt    binary      lcur     ioopt     
       0.1         2         0         1
*DATABASE_JNTFORC
$#      dt    binary      lcur     ioopt     
       0.1         2         0         1
*DATABASE_MATSUM
$#      dt    binary      lcur     ioopt     
       1.0         2         0         1
*DATABASE_NODOUT
$#      dt    binary      lcur     ioopt   option1   option2       
       0.1         2         0         1       0.0         0
*DATABASE_RBDOUT
$#      dt    binary      lcur     ioopt     
       2.0         2         0         1
*DATABASE_RWFORC
$#      dt    binary      lcur     ioopt     
       0.1         2         0         1
*DATABASE_SBTOUT
$#      dt    binary      lcur     ioopt     
       0.1         2         0         1
*DATABASE_SECFORC
$#      dt    binary      lcur     ioopt     
       0.1         2         0         1
*DATABASE_SLEOUT
$#      dt    binary      lcur     ioopt     
       0.1         2         0         1
*DATABASE_SSSTAT
$#      dt    binary      lcur     ioopt     
       0.1         2         0         1
*DATABASE_EXTENT_BINARY
$$   NEIPH     NEIPS    MAXINT    STRFLG    SIGFLG    EPSFLG    RLTFLG    ENGFLG
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
         0         7         3         1         1         1         1         1
$$  CMPFLG    IEVERP    BEAMIP     DCOMP      SHGE     STSSZ    N3THDT   NINTSLD
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         1         0         1         2         2         1         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       0.0         0         0         0                    
$#    dtdt    resplt     neipb     quadr     cubic     
         0         0         0         0         0
$
$
$ =========
$ SET cards
$ =========
$
$
$: Cross-reference summary for SET_PART 6600000
$:---------------------------------------------
$: CONTROL <No label>
$
*SET_PART_LIST_GENERATE
$ This ID 6,600,000 is to be used for psnfail in *CONTROL_SHELL
$#     sid       da1       da2       da3       da4    solver      
   6600000       0.0       0.0       0.0       0.0MECH
$#   b1beg     b1end     b2beg     b2end     b3beg     b3end     b4beg     b4end
         1   6599999   6600001  99999999         0         0         0         0
$
$
*END
