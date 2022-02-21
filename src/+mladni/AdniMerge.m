classdef AdniMerge < handle
    %% e.g.:
    %        OriginalVariableNames                   Var1                           Var2                           Var3                           Var4                           Var5                           Var6                           Var7            
    %     ____________________________    ___________________________    ___________________________    ___________________________    ___________________________    ___________________________    ___________________________    ___________________________
    % 
    %     {'RID'                     }    {[                    230]}    {[                    230]}    {[                    230]}    {[                    230]}    {[                    230]}    {[                    230]}    {[                    230]}
    %     {'PTID'                    }    {'128_S_0230'             }    {'128_S_0230'             }    {'128_S_0230'             }    {'128_S_0230'             }    {'128_S_0230'             }    {'128_S_0230'             }    {'128_S_0230'             }
    %     {'VISCODE'                 }    {'bl'                     }    {'m06'                    }    {'m12'                    }    {'m24'                    }    {'m36'                    }    {'m96'                    }    {'m72'                    }
    %     {'SITE'                    }    {[                    128]}    {[                    128]}    {[                    128]}    {[                    128]}    {[                    128]}    {[                    128]}    {[                    128]}
    %     {'COLPROT'                 }    {'ADNI1'                  }    {'ADNI1'                  }    {'ADNI1'                  }    {'ADNI1'                  }    {'ADNI1'                  }    {'ADNI2'                  }    {'ADNI2'                  }
    %     {'ORIGPROT'                }    {'ADNI1'                  }    {'ADNI1'                  }    {'ADNI1'                  }    {'ADNI1'                  }    {'ADNI1'                  }    {'ADNI1'                  }    {'ADNI1'                  }
    %     {'EXAMDATE'                }    {[2006-03-30             ]}    {[2006-10-18             ]}    {[2007-03-23             ]}    {[2008-04-04             ]}    {[2009-09-16             ]}    {[2014-01-14             ]}    {[2012-01-05             ]}
    %     {'DX_bl'                   }    {'CN'                     }    {'CN'                     }    {'CN'                     }    {'CN'                     }    {'CN'                     }    {'CN'                     }    {'CN'                     }
    %     {'AGE'                     }    {[                     80]}    {[                     80]}    {[                     80]}    {[                     80]}    {[                     80]}    {[                     80]}    {[                     80]}
    %     {'PTGENDER'                }    {'Male'                   }    {'Male'                   }    {'Male'                   }    {'Male'                   }    {'Male'                   }    {'Male'                   }    {'Male'                   }
    %     {'PTEDUCAT'                }    {[                     20]}    {[                     20]}    {[                     20]}    {[                     20]}    {[                     20]}    {[                     20]}    {[                     20]}
    %     {'PTETHCAT'                }    {'Not Hisp/Latino'        }    {'Not Hisp/Latino'        }    {'Not Hisp/Latino'        }    {'Not Hisp/Latino'        }    {'Not Hisp/Latino'        }    {'Not Hisp/Latino'        }    {'Not Hisp/Latino'        }
    %     {'PTRACCAT'                }    {'White'                  }    {'White'                  }    {'White'                  }    {'White'                  }    {'White'                  }    {'White'                  }    {'White'                  }
    %     {'PTMARRY'                 }    {'Married'                }    {'Married'                }    {'Married'                }    {'Married'                }    {'Married'                }    {'Married'                }    {'Married'                }
    %     {'APOE4'                   }    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}
    %     {'FDG'                     }    {[      1.098290000000000]}    {[      1.062490000000000]}    {[      1.057430000000000]}    {[      1.010540000000000]}    {[      1.162290000000000]}    {[      0.785154000000000]}    {[      0.932040000000000]}
    %     {'PIB'                     }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'AV45'                    }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      1.086400000000000]}    {[      1.101900000000000]}
    %     {'ABETA'                   }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'TAU'                     }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'PTAU'                    }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'CDRSB'                   }    {[                      0]}    {[                      1]}    {[      1.500000000000000]}    {[      2.500000000000000]}    {[      5.500000000000000]}    {[                     10]}    {[                      7]}
    %     {'ADAS11'                  }    {[                      8]}    {[     10.330000000000000]}    {[      8.670000000000000]}    {[     13.330000000000000]}    {[                     15]}    {[                     33]}    {[                     17]}
    %     {'ADAS13'                  }    {[                     11]}    {[     13.330000000000000]}    {[     10.670000000000000]}    {[     18.329999999999998]}    {[                     20]}    {[                     46]}    {[                     27]}
    %     {'ADASQ4'                  }    {[                      2]}    {[                      2]}    {[                      2]}    {[                      4]}    {[                      4]}    {[                     10]}    {[                      7]}
    %     {'MMSE'                    }    {[                     29]}    {[                     28]}    {[                     30]}    {[                     27]}    {[                     26]}    {[                     18]}    {[                     22]}
    %     {'RAVLT_immediate'         }    {[                     49]}    {[                     49]}    {[                     54]}    {[                     38]}    {[                     42]}    {[                     23]}    {[                     34]}
    %     {'RAVLT_learning'          }    {[                      5]}    {[                     12]}    {[                      9]}    {[                      8]}    {[                      8]}    {[                      1]}    {[                      4]}
    %     {'RAVLT_forgetting'        }    {[                      1]}    {[                      3]}    {[                      2]}    {[                      2]}    {[                      6]}    {[                      5]}    {[                      5]}
    %     {'RAVLT_perc_forgetting'   }    {[      7.692310000000000]}    {[                     20]}    {[     13.333299999999999]}    {[     18.181799999999999]}    {[                     50]}    {[                    100]}    {[     55.555599999999998]}
    %     {'LDELTOTAL'               }    {[                     10]}    {[                    NaN]}    {[                      8]}    {[                      4]}    {[                      8]}    {[                      0]}    {[                      3]}
    %     {'DIGITSCOR'               }    {[                     38]}    {[                     43]}    {[                     43]}    {[                     36]}    {[                     31]}    {[                    NaN]}    {[                    NaN]}
    %     {'TRABSCOR'                }    {[                    195]}    {[                     74]}    {[                     77]}    {[                    105]}    {[                    128]}    {[                    NaN]}    {[                    283]}
    %     {'FAQ'                     }    {[                      0]}    {[                      2]}    {[                      0]}    {[                      4]}    {[                     21]}    {[                     30]}    {[                     25]}
    %     {'MOCA'                    }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                     13]}    {[                     20]}
    %     {'EcogPtMem'               }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                      3]}    {[                      1]}
    %     {'EcogPtLang'              }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      2.444440000000000]}    {[                      1]}
    %     {'EcogPtVisspat'           }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      2.571430000000000]}    {[      1.142860000000000]}
    %     {'EcogPtPlan'              }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      3.200000000000000]}    {[                      1]}
    %     {'EcogPtOrgan'             }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      3.500000000000000]}    {[                    NaN]}
    %     {'EcogPtDivatt'            }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      3.333330000000000]}    {[                      1]}
    %     {'EcogPtTotal'             }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      2.921050000000000]}    {[      1.029410000000000]}
    %     {'EcogSPMem'               }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                      4]}    {[                      4]}
    %     {'EcogSPLang'              }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      3.888890000000000]}    {[      2.222220000000000]}
    %     {'EcogSPVisspat'           }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      3.857140000000000]}    {[                      3]}
    %     {'EcogSPPlan'              }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                      4]}    {[      3.800000000000000]}
    %     {'EcogSPOrgan'             }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                      4]}    {[                      4]}
    %     {'EcogSPDivatt'            }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                      4]}    {[                    NaN]}
    %     {'EcogSPTotal'             }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[      3.948720000000000]}    {[      3.351350000000000]}
    %     {'FLDSTRENG'               }    {'1.5 Tesla MRI'          }    {'1.5 Tesla MRI'          }    {'1.5 Tesla MRI'          }    {'1.5 Tesla MRI'          }    {0×0 char                 }    {0×0 char                 }    {0×0 char                 }
    %     {'IMAGEUID'                }    {[                 119383]}    {[                 119384]}    {[                  68849]}    {[                 102423]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'Ventricles'              }    {[                  39562]}    {[                  39867]}    {[                  41076]}    {[                  43753]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'Hippocampus'             }    {[                   5700]}    {[                   5888]}    {[                   5787]}    {[                   5830]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'WholeBrain'              }    {[                1051050]}    {[                1076940]}    {[                1064480]}    {[                1052260]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'Entorhinal'              }    {[                   2366]}    {[                   2295]}    {[                   2510]}    {[                   2735]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'Fusiform'                }    {[                  16394]}    {[                  16664]}    {[                  16898]}    {[                  16465]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'MidTemp'                 }    {[                  22334]}    {[                  22834]}    {[                  21972]}    {[                  22467]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'ICV'                     }    {[                1714030]}    {[                1739610]}    {[                1742490]}    {[                1745920]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'DX'                      }    {'CN'                     }    {'CN'                     }    {'CN'                     }    {'CN'                     }    {'Dementia'               }    {'Dementia'               }    {'Dementia'               }
    %     {'mPACCdigit'              }    {[     -1.392850000000000]}    {[     -1.084710000000000]}    {[     -0.602291000000000]}    {[     -6.295530000000000]}    {[     -6.487120000000000]}    {[    -23.953900000000001]}    {[    -15.737299999999999]}
    %     {'mPACCtrailsB'            }    {[     -2.901980000000000]}    {[     -0.683706000000000]}    {[     -0.395445000000000]}    {[     -6.136900000000000]}    {[     -6.308540000000000]}    {[    -23.953900000000001]}    {[    -14.959400000000000]}
    %     {'EXAMDATE_bl'             }    {[2006-03-30             ]}    {[2006-03-30             ]}    {[2006-03-30             ]}    {[2006-03-30             ]}    {[2006-03-30             ]}    {[2006-03-30             ]}    {[2006-03-30             ]}
    %     {'CDRSB_bl'                }    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}
    %     {'ADAS11_bl'               }    {[                      8]}    {[                      8]}    {[                      8]}    {[                      8]}    {[                      8]}    {[                      8]}    {[                      8]}
    %     {'ADAS13_bl'               }    {[                     11]}    {[                     11]}    {[                     11]}    {[                     11]}    {[                     11]}    {[                     11]}    {[                     11]}
    %     {'ADASQ4_bl'               }    {[                      2]}    {[                      2]}    {[                      2]}    {[                      2]}    {[                      2]}    {[                      2]}    {[                      2]}
    %     {'MMSE_bl'                 }    {[                     29]}    {[                     29]}    {[                     29]}    {[                     29]}    {[                     29]}    {[                     29]}    {[                     29]}
    %     {'RAVLT_immediate_bl'      }    {[                     49]}    {[                     49]}    {[                     49]}    {[                     49]}    {[                     49]}    {[                     49]}    {[                     49]}
    %     {'RAVLT_learning_bl'       }    {[                      5]}    {[                      5]}    {[                      5]}    {[                      5]}    {[                      5]}    {[                      5]}    {[                      5]}
    %     {'RAVLT_forgetting_bl'     }    {[                      1]}    {[                      1]}    {[                      1]}    {[                      1]}    {[                      1]}    {[                      1]}    {[                      1]}
    %     {'RAVLT_perc_forgetting_bl'}    {[      7.692310000000000]}    {[      7.692310000000000]}    {[      7.692310000000000]}    {[      7.692310000000000]}    {[      7.692310000000000]}    {[      7.692310000000000]}    {[      7.692310000000000]}
    %     {'LDELTOTAL_BL'            }    {[                     10]}    {[                     10]}    {[                     10]}    {[                     10]}    {[                     10]}    {[                     10]}    {[                     10]}
    %     {'DIGITSCOR_bl'            }    {[                     38]}    {[                     38]}    {[                     38]}    {[                     38]}    {[                     38]}    {[                     38]}    {[                     38]}
    %     {'TRABSCOR_bl'             }    {[                    195]}    {[                    195]}    {[                    195]}    {[                    195]}    {[                    195]}    {[                    195]}    {[                    195]}
    %     {'FAQ_bl'                  }    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}    {[                      0]}
    %     {'mPACCdigit_bl'           }    {[     -1.392850000000000]}    {[     -1.392850000000000]}    {[     -1.392850000000000]}    {[     -1.392850000000000]}    {[     -1.392850000000000]}    {[     -1.392850000000000]}    {[     -1.392850000000000]}
    %     {'mPACCtrailsB_bl'         }    {[     -2.901980000000000]}    {[     -2.901980000000000]}    {[     -2.901980000000000]}    {[     -2.901980000000000]}    {[     -2.901980000000000]}    {[     -2.901980000000000]}    {[     -2.901980000000000]}
    %     {'FLDSTRENG_bl'            }    {'1.5 Tesla MRI'          }    {'1.5 Tesla MRI'          }    {'1.5 Tesla MRI'          }    {'1.5 Tesla MRI'          }    {'1.5 Tesla MRI'          }    {'1.5 Tesla MRI'          }    {'1.5 Tesla MRI'          }
    %     {'Ventricles_bl'           }    {[                  39562]}    {[                  39562]}    {[                  39562]}    {[                  39562]}    {[                  39562]}    {[                  39562]}    {[                  39562]}
    %     {'Hippocampus_bl'          }    {[                   5700]}    {[                   5700]}    {[                   5700]}    {[                   5700]}    {[                   5700]}    {[                   5700]}    {[                   5700]}
    %     {'WholeBrain_bl'           }    {[                1051050]}    {[                1051050]}    {[                1051050]}    {[                1051050]}    {[                1051050]}    {[                1051050]}    {[                1051050]}
    %     {'Entorhinal_bl'           }    {[                   2366]}    {[                   2366]}    {[                   2366]}    {[                   2366]}    {[                   2366]}    {[                   2366]}    {[                   2366]}
    %     {'Fusiform_bl'             }    {[                  16394]}    {[                  16394]}    {[                  16394]}    {[                  16394]}    {[                  16394]}    {[                  16394]}    {[                  16394]}
    %     {'MidTemp_bl'              }    {[                  22334]}    {[                  22334]}    {[                  22334]}    {[                  22334]}    {[                  22334]}    {[                  22334]}    {[                  22334]}
    %     {'ICV_bl'                  }    {[                1714030]}    {[                1714030]}    {[                1714030]}    {[                1714030]}    {[                1714030]}    {[                1714030]}    {[                1714030]}
    %     {'MOCA_bl'                 }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogPtMem_bl'            }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogPtLang_bl'           }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogPtVisspat_bl'        }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogPtPlan_bl'           }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogPtOrgan_bl'          }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogPtDivatt_bl'         }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogPtTotal_bl'          }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogSPMem_bl'            }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogSPLang_bl'           }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogSPVisspat_bl'        }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogSPPlan_bl'           }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogSPOrgan_bl'          }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogSPDivatt_bl'         }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'EcogSPTotal_bl'          }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'ABETA_bl'                }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'TAU_bl'                  }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'PTAU_bl'                 }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'FDG_bl'                  }    {[      1.098290000000000]}    {[      1.098290000000000]}    {[      1.098290000000000]}    {[      1.098290000000000]}    {[      1.098290000000000]}    {[      1.098290000000000]}    {[      1.098290000000000]}
    %     {'PIB_bl'                  }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'AV45_bl'                 }    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}    {[                    NaN]}
    %     {'Years_bl'                }    {[                      0]}    {[      0.553046000000000]}    {[      0.980151000000000]}    {[      2.015060000000000]}    {[      3.466120000000000]}    {[      7.794660000000000]}    {[      5.768650000000000]}
    %     {'Month_bl'                }    {[                      0]}    {[      6.622950000000000]}    {[     11.737700000000000]}    {[     24.131100000000000]}    {[     41.508200000000002]}    {[     93.344300000000004]}    {[     69.081999999999994]}
    %     {'Month'                   }    {[                      0]}    {[                      6]}    {[                     12]}    {[                     24]}    {[                     42]}    {[                     96]}    {[                     72]}
    %     {'M'                       }    {[                      0]}    {[                      6]}    {[                     12]}    {[                     24]}    {[                     36]}    {[                     96]}    {[                     72]}
    %     {'update_stamp'            }    {[2021-07-31 06:16:38.000]}    {[2021-07-31 06:16:38.000]}    {[2021-07-31 06:16:38.000]}    {[2021-07-31 06:16:38.000]}    {[2021-07-31 06:16:38.000]}    {[2021-07-31 06:16:38.000]}    {[2021-07-31 06:16:38.000]}
    %     
    %  Created 26-Dec-2021 19:07:31 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1809720 (R2021b) Update 1 for MACI64.  Copyright 2021 John J. Lee.
    
    properties
        home
    end

    properties (Constant)
        categories = {'CN', 'EMCI', 'MCI', 'LMCI', 'SMC', 'AD'};
    end

    properties (Dependent)
        merge_file
        dict_file
        subjects
    end

    methods

        %% GET

        function g = get.merge_file(this)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "ADNIMERGE.csv");
        end
        function g = get.dict_file(this)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "ADNIMERGE_DICT.csv");
        end
        function g = get.subjects(this)
            g = this.subjects_;
        end       

        %%

        function d = dict(this, varargin)
            d = this.dict_;
        end
        function t = table(this, varargin)
            t = this.merge_;
        end
        function t = table_dCDRSB(this, varargin)
            if ~isempty(this.dCDRSB_)
                t = this.dCDRSB_;
                return
            end

            % lazy init
            s = this.subjects_(1);
            t_ = this.merge_;
            t__ = t_(t_.RID == s & ~isnan(t_.FDG), :); % pick subject subtable
            d = max(t__.CDRSB) - min(t__.CDRSB_bl); % scalar dCDRSB
            t = table(s, d, 'VariableNames', {'RID', 'dCDRSB'});
            for ui = 2:length(this.subjects_)
                s = this.subjects_(ui);
                u__ = t_(t_.RID == s & ~isnan(t_.FDG), :); 
                d = max(u__.CDRSB) - min(u__.CDRSB_bl);
                try
                    t = [t; table(s, d, 'VariableNames', {'RID', 'dCDRSB'})];
                catch
                    %fprintf("no reasonable subtable for RID ~ %g\n", this.subjects_(ui))
                end
            end
            this.dCDRSB_ = t;
        end
        function t = table_CDRSB_from_bl(this, varargin)
            %  bl_thresh ~ 0 => 411 subjects
            %              1 => 792
            %              2 => 1074
            %              3 => 1204
            %              4 => 1277
            %              5 => 1360
            %  total avail. subjects ~ 1855

            ip = inputParser;
            addOptional(ip, 'bl_thresh', 1, @isscalar) % select subjects with baseline <= bl_thresh
            addParameter(ip, 'use_cached', false, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;

            if ipr.use_cached && ~isempty(this.CDRSB_bl_)
                t = this.CDRSB_bl_;
                return
            end

            % lazy init
            s = this.subjects_(1);
            t_ = this.merge_;
            t_ = t_(t_.CDRSB_bl <= ipr.bl_thresh, :);
            t__ = t_(t_.RID == s & ~isnan(t_.FDG), :); % pick subject subtable
            d = max(t__.CDRSB); % scalar dCDRSB
            t = table(s, d, 'VariableNames', {'RID', 'CDRSB'});
            for ui = 2:length(this.subjects_)
                s = this.subjects_(ui);
                u__ = t_(t_.RID == s & ~isnan(t_.FDG), :); 
                d = max(u__.CDRSB);
                try
                    t = [t; table(s, d, 'VariableNames', {'RID', 'CDRSB'})];
                catch
                    %fprintf("no reasonable subtable for RID ~ %g\n", this.subjects_(ui))
                end
            end
            this.CDRSB_bl_ = t;
        end
        function v = variableNames(this)
            v = this.merge_.Properties.VariableNames;
        end

        function this = AdniMerge(varargin)
            ip = inputParser;
            addParameter(ip, "home", ...
                fullfile(getenv("SINGULARITY_HOME"), "ADNI", "testing", ""), @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.home = ipr.home;
            
            cd(this.home)
            this.merge_ = readtable(this.merge_file);
            dict_full_ = readtable(this.dict_file);
            this.dict_ = table(dict_full_.FLDNAME, dict_full_.TEXT, 'VariableNames', {'field', 'description'});
            this.subjects_ = unique(this.merge_.RID);
            
        end
    end

    %% PROTECTED

    properties (Access = protected)        
        CDRSB_bl_
        dCDRSB_
        dict_
        lastscan_
        merge_
        subjects_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
