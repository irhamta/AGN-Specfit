FUNCTION gfit_eval_D0, x, par, EXPR=_expr
  COMPILE_OPT IDL2
  ON_ERROR, !glib.on_error
  COMMON gfit
  COMMON gfit_private
  currentDataSet = 0

  CONTINUUM = qsfit_comp_sbpowerlaw(x, par[0], par[1], par[2], par[3], par[4])
  GALAXY = qsfit_comp_galaxytemplate(x, par[5], _extra=gfit.comp.GALAXY.opt)
  BALMER = qsfit_comp_balmer(x, par[6], par[7], par[8], par[9], par[10], par[11])
  IRONUV = qsfit_comp_ironuv(x, par[12], par[13])
  IRONOPT = qsfit_comp_ironoptical(x, par[14], par[15], par[16], par[17])
  BR_SIIV_1400 = 0
  BR_CIV_1549 = 0
  BR_CIII_1909 = 0
  BR_MGII_2798 = qsfit_comp_emline(x, par[30], par[31], par[32], par[33])
  NA_NEVI_3426 = qsfit_comp_emline(x, par[34], par[35], par[36], par[37])
  NA_OII_3727 = qsfit_comp_emline(x, par[38], par[39], par[40], par[41])
  NA_NEIII_3869 = qsfit_comp_emline(x, par[42], par[43], par[44], par[45])
  BR_HD = qsfit_comp_emline(x, par[46], par[47], par[48], par[49])
  BR_HG = qsfit_comp_emline(x, par[50], par[51], par[52], par[53])
  NA_HB = qsfit_comp_emline(x, par[54], par[55], par[56], par[57])
  BR_HB = qsfit_comp_emline(x, par[58], par[59], par[60], par[61])
  NA_OIII_4959 = qsfit_comp_emline(x, par[62], par[63], par[64], par[65])
  NA_OIII_5007 = qsfit_comp_emline(x, par[66], par[67], par[68], par[69])
  BR_HEI_5876 = 0
  NA_NII_6549 = 0
  NA_HA = 0
  BR_HA = 0
  NA_NII_6583 = 0
  NA_SII_6716 = 0
  NA_SII_6731 = 0
  UNK1 = qsfit_comp_emline(x, par[98], par[99], par[100], par[101])
  UNK2 = qsfit_comp_emline(x, par[102], par[103], par[104], par[105])
  UNK3 = qsfit_comp_emline(x, par[106], par[107], par[108], par[109])
  UNK4 = qsfit_comp_emline(x, par[110], par[111], par[112], par[113])
  UNK5 = qsfit_comp_emline(x, par[114], par[115], par[116], par[117])
  UNK6 = qsfit_comp_emline(x, par[118], par[119], par[120], par[121])
  UNK7 = qsfit_comp_emline(x, par[122], par[123], par[124], par[125])
  UNK8 = qsfit_comp_emline(x, par[126], par[127], par[128], par[129])
  UNK9 = qsfit_comp_emline(x, par[130], par[131], par[132], par[133])
  UNK10 = qsfit_comp_emline(x, par[134], par[135], par[136], par[137])
  LINE_HA_BASE = 0

  MODEL = CONTINUUM + GALAXY + BALMER + IRONUV + IRONOPT + BR_SIIV_1400 + BR_CIV_1549 + BR_CIII_1909 + BR_MGII_2798 + NA_NEVI_3426 + NA_OII_3727 + NA_NEIII_3869 + BR_HD + BR_HG + NA_HB + BR_HB + NA_OIII_4959 + NA_OIII_5007 + BR_HEI_5876 + NA_NII_6549 + NA_HA + BR_HA + NA_NII_6583 + NA_SII_6716 + NA_SII_6731 + UNK1 + UNK2 + UNK3 + UNK4 + UNK5 + UNK6 + UNK7 + UNK8 + UNK9 + UNK10 + LINE_HA_BASE

  IF (flag_evalAllExpr) THEN BEGIN
    EXPR_CONTINUUM = continuum
    EXPR_GALAXY = galaxy
    EXPR_CONTGALAXY = continuum + galaxy
    EXPR_BALMER = balmer
    EXPR_IRON = ironuv + ironopt
    EXPR_BROADLINES = line_Ha_base + br_SiIV_1400 + br_CIV_1549 + br_CIII_1909 + br_MgII_2798 + br_Hd + br_Hg + br_Hb + br_HeI_5876 + br_Ha
    EXPR_NARROWLINES = na_NeVI_3426 + na_OII_3727 + na_NeIII_3869 + na_Hb + na_OIII_4959 + na_OIII_5007 + na_NII_6549 + na_Ha + na_NII_6583 + na_SII_6716 + na_SII_6731
    EXPR_UNKNOWN = unk1 + unk2 + unk3 + unk4 + unk5 + unk6 + unk7 + unk8 + unk9 + unk10
    _expr = { MODEL: MODEL, EXPR_CONTINUUM: EXPR_CONTINUUM, EXPR_GALAXY: EXPR_GALAXY, EXPR_CONTGALAXY: EXPR_CONTGALAXY, EXPR_BALMER: EXPR_BALMER, EXPR_IRON: EXPR_IRON, EXPR_BROADLINES: EXPR_BROADLINES, EXPR_NARROWLINES: EXPR_NARROWLINES, EXPR_UNKNOWN: EXPR_UNKNOWN}
  ENDIF

  RETURN, MODEL
END



FUNCTION mpfit_eval_model, par, EXPR=expr
  COMPILE_OPT IDL2
  ON_ERROR, !glib.on_error
  COMMON gfit
  COMMON gfit_private

  flag_evalAllExpr = ARG_PRESENT(expr)
  gfit.cmp.D0.m = gfit_eval_D0(gfit.data.D0.x[dataNoticed.(0)], par, EXPR=expr_D0)

  IF (flag_evalAllExpr) THEN BEGIN
    expr = { D0: expr_D0}
  ENDIF

  wdev_D0 = gfit_weighted_dev_gauss( gfit.cmp.D0 )
  RETURN, [wdev_D0]
END



PRO mpfit_eval_model0
END

