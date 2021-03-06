/*Header-MicMac-eLiSe-25/06/2007

    MicMac : Multi Image Correspondances par Methodes Automatiques de Correlation
    eLiSe  : ELements of an Image Software Environnement

    www.micmac.ign.fr

   
    Copyright : Institut Geographique National
    Author : Marc Pierrot Deseilligny
    Contributors : Gregoire Maillet, Didier Boldo.

[1] M. Pierrot-Deseilligny, N. Paparoditis.
    "A multiresolution and optimization-based image matching approach:
    An application to surface reconstruction from SPOT5-HRS stereo imagery."
    In IAPRS vol XXXVI-1/W41 in ISPRS Workshop On Topographic Mapping From Space
    (With Special Emphasis on Small Satellites), Ankara, Turquie, 02-2006.

[2] M. Pierrot-Deseilligny, "MicMac, un lociel de mise en correspondance
    d'images, adapte au contexte geograhique" to appears in 
    Bulletin d'information de l'Institut Geographique National, 2007.

Francais :

   MicMac est un logiciel de mise en correspondance d'image adapte 
   au contexte de recherche en information geographique. Il s'appuie sur
   la bibliotheque de manipulation d'image eLiSe. Il est distibue sous la
   licences Cecill-B.  Voir en bas de fichier et  http://www.cecill.info.


English :

    MicMac is an open source software specialized in image matching
    for research in geographic information. MicMac is built on the
    eLiSe image library. MicMac is governed by the  "Cecill-B licence".
    See below and http://www.cecill.info.

Header-MicMac-eLiSe-25/06/2007*/
// File Automatically generated by eLiSe

#include "general/all.h"
#include "private/all.h"
#include "cEqCorrelGrid_9_Im2Var.h"


cEqCorrelGrid_9_Im2Var::cEqCorrelGrid_9_Im2Var():
    cElCompiledFonc(9)
{
   AddIntRef (cIncIntervale("ZCorrel0",0,1));
   AddIntRef (cIncIntervale("ZCorrel1",1,2));
   AddIntRef (cIncIntervale("ZCorrel2",2,3));
   AddIntRef (cIncIntervale("ZCorrel3",3,4));
   AddIntRef (cIncIntervale("ZCorrel4",4,5));
   AddIntRef (cIncIntervale("ZCorrel5",5,6));
   AddIntRef (cIncIntervale("ZCorrel6",6,7));
   AddIntRef (cIncIntervale("ZCorrel7",7,8));
   AddIntRef (cIncIntervale("ZCorrel8",8,9));
   Close(false);
}



void cEqCorrelGrid_9_Im2Var::ComputeVal()
{
   double tmp0_ = mLocGr1_0+mLocGr1_1;
   double tmp1_ = tmp0_+mLocGr1_2;
   double tmp2_ = tmp1_+mLocGr1_3;
   double tmp3_ = tmp2_+mLocGr1_4;
   double tmp4_ = tmp3_+mLocGr1_5;
   double tmp5_ = tmp4_+mLocGr1_6;
   double tmp6_ = tmp5_+mLocGr1_7;
   double tmp7_ = tmp6_+mLocGr1_8;
   double tmp8_ = (tmp7_)/9;
   double tmp9_ = mCompCoord[0];
   double tmp10_ = mLocDGr2Dz0*tmp9_;
   double tmp11_ = mLocGr2of0_0+tmp10_;
   double tmp12_ = mCompCoord[1];
   double tmp13_ = mLocDGr2Dz1*tmp12_;
   double tmp14_ = mLocGr2of0_1+tmp13_;
   double tmp15_ = mCompCoord[2];
   double tmp16_ = mLocDGr2Dz2*tmp15_;
   double tmp17_ = mLocGr2of0_2+tmp16_;
   double tmp18_ = mCompCoord[3];
   double tmp19_ = mLocDGr2Dz3*tmp18_;
   double tmp20_ = mLocGr2of0_3+tmp19_;
   double tmp21_ = mCompCoord[4];
   double tmp22_ = mLocDGr2Dz4*tmp21_;
   double tmp23_ = mLocGr2of0_4+tmp22_;
   double tmp24_ = mCompCoord[5];
   double tmp25_ = mLocDGr2Dz5*tmp24_;
   double tmp26_ = mLocGr2of0_5+tmp25_;
   double tmp27_ = mCompCoord[6];
   double tmp28_ = mLocDGr2Dz6*tmp27_;
   double tmp29_ = mLocGr2of0_6+tmp28_;
   double tmp30_ = mCompCoord[7];
   double tmp31_ = mLocDGr2Dz7*tmp30_;
   double tmp32_ = mLocGr2of0_7+tmp31_;
   double tmp33_ = mCompCoord[8];
   double tmp34_ = mLocDGr2Dz8*tmp33_;
   double tmp35_ = mLocGr2of0_8+tmp34_;
   double tmp36_ = tmp11_+tmp14_;
   double tmp37_ = tmp36_+tmp17_;
   double tmp38_ = tmp37_+tmp20_;
   double tmp39_ = tmp38_+tmp23_;
   double tmp40_ = tmp39_+tmp26_;
   double tmp41_ = tmp40_+tmp29_;
   double tmp42_ = tmp41_+tmp32_;
   double tmp43_ = tmp42_+tmp35_;
   double tmp44_ = (tmp43_)/9;
   double tmp45_ = ElSquare(mLocGr1_0);
   double tmp46_ = ElSquare(mLocGr1_1);
   double tmp47_ = tmp45_+tmp46_;
   double tmp48_ = ElSquare(mLocGr1_2);
   double tmp49_ = tmp47_+tmp48_;
   double tmp50_ = ElSquare(mLocGr1_3);
   double tmp51_ = tmp49_+tmp50_;
   double tmp52_ = ElSquare(mLocGr1_4);
   double tmp53_ = tmp51_+tmp52_;
   double tmp54_ = ElSquare(mLocGr1_5);
   double tmp55_ = tmp53_+tmp54_;
   double tmp56_ = ElSquare(mLocGr1_6);
   double tmp57_ = tmp55_+tmp56_;
   double tmp58_ = ElSquare(mLocGr1_7);
   double tmp59_ = tmp57_+tmp58_;
   double tmp60_ = ElSquare(mLocGr1_8);
   double tmp61_ = tmp59_+tmp60_;
   double tmp62_ = (tmp61_)/9;
   double tmp63_ = ElSquare(tmp8_);
   double tmp64_ = tmp62_-tmp63_;
   double tmp65_ = tmp64_+0.001000;
   double tmp66_ = sqrt(tmp65_);
   double tmp67_ = ElSquare(tmp11_);
   double tmp68_ = ElSquare(tmp14_);
   double tmp69_ = tmp67_+tmp68_;
   double tmp70_ = ElSquare(tmp17_);
   double tmp71_ = tmp69_+tmp70_;
   double tmp72_ = ElSquare(tmp20_);
   double tmp73_ = tmp71_+tmp72_;
   double tmp74_ = ElSquare(tmp23_);
   double tmp75_ = tmp73_+tmp74_;
   double tmp76_ = ElSquare(tmp26_);
   double tmp77_ = tmp75_+tmp76_;
   double tmp78_ = ElSquare(tmp29_);
   double tmp79_ = tmp77_+tmp78_;
   double tmp80_ = ElSquare(tmp32_);
   double tmp81_ = tmp79_+tmp80_;
   double tmp82_ = ElSquare(tmp35_);
   double tmp83_ = tmp81_+tmp82_;
   double tmp84_ = (tmp83_)/9;
   double tmp85_ = ElSquare(tmp44_);
   double tmp86_ = tmp84_-tmp85_;
   double tmp87_ = tmp86_+0.001000;
   double tmp88_ = sqrt(tmp87_);

  mVal[0] = (mLocGr1_0-tmp8_)/tmp66_-((tmp11_)-tmp44_)/tmp88_;

  mVal[1] = (mLocGr1_1-tmp8_)/tmp66_-((tmp14_)-tmp44_)/tmp88_;

  mVal[2] = (mLocGr1_2-tmp8_)/tmp66_-((tmp17_)-tmp44_)/tmp88_;

  mVal[3] = (mLocGr1_3-tmp8_)/tmp66_-((tmp20_)-tmp44_)/tmp88_;

  mVal[4] = (mLocGr1_4-tmp8_)/tmp66_-((tmp23_)-tmp44_)/tmp88_;

  mVal[5] = (mLocGr1_5-tmp8_)/tmp66_-((tmp26_)-tmp44_)/tmp88_;

  mVal[6] = (mLocGr1_6-tmp8_)/tmp66_-((tmp29_)-tmp44_)/tmp88_;

  mVal[7] = (mLocGr1_7-tmp8_)/tmp66_-((tmp32_)-tmp44_)/tmp88_;

  mVal[8] = (mLocGr1_8-tmp8_)/tmp66_-((tmp35_)-tmp44_)/tmp88_;

}


void cEqCorrelGrid_9_Im2Var::ComputeValDeriv()
{
   double tmp0_ = mLocGr1_0+mLocGr1_1;
   double tmp1_ = tmp0_+mLocGr1_2;
   double tmp2_ = tmp1_+mLocGr1_3;
   double tmp3_ = tmp2_+mLocGr1_4;
   double tmp4_ = tmp3_+mLocGr1_5;
   double tmp5_ = tmp4_+mLocGr1_6;
   double tmp6_ = tmp5_+mLocGr1_7;
   double tmp7_ = tmp6_+mLocGr1_8;
   double tmp8_ = (tmp7_)/9;
   double tmp9_ = mCompCoord[0];
   double tmp10_ = mLocDGr2Dz0*tmp9_;
   double tmp11_ = mLocGr2of0_0+tmp10_;
   double tmp12_ = mCompCoord[1];
   double tmp13_ = mLocDGr2Dz1*tmp12_;
   double tmp14_ = mLocGr2of0_1+tmp13_;
   double tmp15_ = mCompCoord[2];
   double tmp16_ = mLocDGr2Dz2*tmp15_;
   double tmp17_ = mLocGr2of0_2+tmp16_;
   double tmp18_ = mCompCoord[3];
   double tmp19_ = mLocDGr2Dz3*tmp18_;
   double tmp20_ = mLocGr2of0_3+tmp19_;
   double tmp21_ = mCompCoord[4];
   double tmp22_ = mLocDGr2Dz4*tmp21_;
   double tmp23_ = mLocGr2of0_4+tmp22_;
   double tmp24_ = mCompCoord[5];
   double tmp25_ = mLocDGr2Dz5*tmp24_;
   double tmp26_ = mLocGr2of0_5+tmp25_;
   double tmp27_ = mCompCoord[6];
   double tmp28_ = mLocDGr2Dz6*tmp27_;
   double tmp29_ = mLocGr2of0_6+tmp28_;
   double tmp30_ = mCompCoord[7];
   double tmp31_ = mLocDGr2Dz7*tmp30_;
   double tmp32_ = mLocGr2of0_7+tmp31_;
   double tmp33_ = mCompCoord[8];
   double tmp34_ = mLocDGr2Dz8*tmp33_;
   double tmp35_ = mLocGr2of0_8+tmp34_;
   double tmp36_ = tmp11_+tmp14_;
   double tmp37_ = tmp36_+tmp17_;
   double tmp38_ = tmp37_+tmp20_;
   double tmp39_ = tmp38_+tmp23_;
   double tmp40_ = tmp39_+tmp26_;
   double tmp41_ = tmp40_+tmp29_;
   double tmp42_ = tmp41_+tmp32_;
   double tmp43_ = tmp42_+tmp35_;
   double tmp44_ = (tmp43_)/9;
   double tmp45_ = ElSquare(tmp11_);
   double tmp46_ = ElSquare(tmp14_);
   double tmp47_ = tmp45_+tmp46_;
   double tmp48_ = ElSquare(tmp17_);
   double tmp49_ = tmp47_+tmp48_;
   double tmp50_ = ElSquare(tmp20_);
   double tmp51_ = tmp49_+tmp50_;
   double tmp52_ = ElSquare(tmp23_);
   double tmp53_ = tmp51_+tmp52_;
   double tmp54_ = ElSquare(tmp26_);
   double tmp55_ = tmp53_+tmp54_;
   double tmp56_ = ElSquare(tmp29_);
   double tmp57_ = tmp55_+tmp56_;
   double tmp58_ = ElSquare(tmp32_);
   double tmp59_ = tmp57_+tmp58_;
   double tmp60_ = ElSquare(tmp35_);
   double tmp61_ = tmp59_+tmp60_;
   double tmp62_ = (tmp61_)/9;
   double tmp63_ = ElSquare(tmp44_);
   double tmp64_ = tmp62_-tmp63_;
   double tmp65_ = tmp64_+0.001000;
   double tmp66_ = sqrt(tmp65_);
   double tmp67_ = (tmp11_)-tmp44_;
   double tmp68_ = ElSquare(9);
   double tmp69_ = mLocDGr2Dz0*9;
   double tmp70_ = (tmp69_)/tmp68_;
   double tmp71_ = mLocDGr2Dz1*9;
   double tmp72_ = (tmp71_)/tmp68_;
   double tmp73_ = ElSquare(tmp66_);
   double tmp74_ = mLocDGr2Dz2*9;
   double tmp75_ = (tmp74_)/tmp68_;
   double tmp76_ = mLocDGr2Dz3*9;
   double tmp77_ = (tmp76_)/tmp68_;
   double tmp78_ = mLocDGr2Dz4*9;
   double tmp79_ = (tmp78_)/tmp68_;
   double tmp80_ = mLocDGr2Dz5*9;
   double tmp81_ = (tmp80_)/tmp68_;
   double tmp82_ = mLocDGr2Dz6*9;
   double tmp83_ = (tmp82_)/tmp68_;
   double tmp84_ = mLocDGr2Dz7*9;
   double tmp85_ = (tmp84_)/tmp68_;
   double tmp86_ = mLocDGr2Dz8*9;
   double tmp87_ = (tmp86_)/tmp68_;
   double tmp88_ = ElSquare(mLocGr1_0);
   double tmp89_ = ElSquare(mLocGr1_1);
   double tmp90_ = tmp88_+tmp89_;
   double tmp91_ = ElSquare(mLocGr1_2);
   double tmp92_ = tmp90_+tmp91_;
   double tmp93_ = ElSquare(mLocGr1_3);
   double tmp94_ = tmp92_+tmp93_;
   double tmp95_ = ElSquare(mLocGr1_4);
   double tmp96_ = tmp94_+tmp95_;
   double tmp97_ = ElSquare(mLocGr1_5);
   double tmp98_ = tmp96_+tmp97_;
   double tmp99_ = ElSquare(mLocGr1_6);
   double tmp100_ = tmp98_+tmp99_;
   double tmp101_ = ElSquare(mLocGr1_7);
   double tmp102_ = tmp100_+tmp101_;
   double tmp103_ = ElSquare(mLocGr1_8);
   double tmp104_ = tmp102_+tmp103_;
   double tmp105_ = (tmp104_)/9;
   double tmp106_ = ElSquare(tmp8_);
   double tmp107_ = tmp105_-tmp106_;
   double tmp108_ = tmp107_+0.001000;
   double tmp109_ = sqrt(tmp108_);
   double tmp110_ = (tmp14_)-tmp44_;
   double tmp111_ = 2*mLocDGr2Dz0;
   double tmp112_ = tmp111_*(tmp11_);
   double tmp113_ = tmp112_*9;
   double tmp114_ = (tmp113_)/tmp68_;
   double tmp115_ = 2*(tmp70_);
   double tmp116_ = tmp115_*(tmp44_);
   double tmp117_ = tmp114_-tmp116_;
   double tmp118_ = 0.500000*(tmp117_);
   double tmp119_ = (tmp118_)/tmp66_;
   double tmp120_ = 2*mLocDGr2Dz1;
   double tmp121_ = tmp120_*(tmp14_);
   double tmp122_ = tmp121_*9;
   double tmp123_ = (tmp122_)/tmp68_;
   double tmp124_ = 2*(tmp72_);
   double tmp125_ = tmp124_*(tmp44_);
   double tmp126_ = tmp123_-tmp125_;
   double tmp127_ = 0.500000*(tmp126_);
   double tmp128_ = (tmp127_)/tmp66_;
   double tmp129_ = -(tmp75_);
   double tmp130_ = tmp129_*tmp66_;
   double tmp131_ = 2*mLocDGr2Dz2;
   double tmp132_ = tmp131_*(tmp17_);
   double tmp133_ = tmp132_*9;
   double tmp134_ = (tmp133_)/tmp68_;
   double tmp135_ = 2*(tmp75_);
   double tmp136_ = tmp135_*(tmp44_);
   double tmp137_ = tmp134_-tmp136_;
   double tmp138_ = 0.500000*(tmp137_);
   double tmp139_ = (tmp138_)/tmp66_;
   double tmp140_ = -(tmp77_);
   double tmp141_ = tmp140_*tmp66_;
   double tmp142_ = 2*mLocDGr2Dz3;
   double tmp143_ = tmp142_*(tmp20_);
   double tmp144_ = tmp143_*9;
   double tmp145_ = (tmp144_)/tmp68_;
   double tmp146_ = 2*(tmp77_);
   double tmp147_ = tmp146_*(tmp44_);
   double tmp148_ = tmp145_-tmp147_;
   double tmp149_ = 0.500000*(tmp148_);
   double tmp150_ = (tmp149_)/tmp66_;
   double tmp151_ = -(tmp79_);
   double tmp152_ = tmp151_*tmp66_;
   double tmp153_ = 2*mLocDGr2Dz4;
   double tmp154_ = tmp153_*(tmp23_);
   double tmp155_ = tmp154_*9;
   double tmp156_ = (tmp155_)/tmp68_;
   double tmp157_ = 2*(tmp79_);
   double tmp158_ = tmp157_*(tmp44_);
   double tmp159_ = tmp156_-tmp158_;
   double tmp160_ = 0.500000*(tmp159_);
   double tmp161_ = (tmp160_)/tmp66_;
   double tmp162_ = -(tmp81_);
   double tmp163_ = tmp162_*tmp66_;
   double tmp164_ = 2*mLocDGr2Dz5;
   double tmp165_ = tmp164_*(tmp26_);
   double tmp166_ = tmp165_*9;
   double tmp167_ = (tmp166_)/tmp68_;
   double tmp168_ = 2*(tmp81_);
   double tmp169_ = tmp168_*(tmp44_);
   double tmp170_ = tmp167_-tmp169_;
   double tmp171_ = 0.500000*(tmp170_);
   double tmp172_ = (tmp171_)/tmp66_;
   double tmp173_ = -(tmp83_);
   double tmp174_ = tmp173_*tmp66_;
   double tmp175_ = 2*mLocDGr2Dz6;
   double tmp176_ = tmp175_*(tmp29_);
   double tmp177_ = tmp176_*9;
   double tmp178_ = (tmp177_)/tmp68_;
   double tmp179_ = 2*(tmp83_);
   double tmp180_ = tmp179_*(tmp44_);
   double tmp181_ = tmp178_-tmp180_;
   double tmp182_ = 0.500000*(tmp181_);
   double tmp183_ = (tmp182_)/tmp66_;
   double tmp184_ = -(tmp85_);
   double tmp185_ = tmp184_*tmp66_;
   double tmp186_ = 2*mLocDGr2Dz7;
   double tmp187_ = tmp186_*(tmp32_);
   double tmp188_ = tmp187_*9;
   double tmp189_ = (tmp188_)/tmp68_;
   double tmp190_ = 2*(tmp85_);
   double tmp191_ = tmp190_*(tmp44_);
   double tmp192_ = tmp189_-tmp191_;
   double tmp193_ = 0.500000*(tmp192_);
   double tmp194_ = (tmp193_)/tmp66_;
   double tmp195_ = -(tmp87_);
   double tmp196_ = tmp195_*tmp66_;
   double tmp197_ = 2*mLocDGr2Dz8;
   double tmp198_ = tmp197_*(tmp35_);
   double tmp199_ = tmp198_*9;
   double tmp200_ = (tmp199_)/tmp68_;
   double tmp201_ = 2*(tmp87_);
   double tmp202_ = tmp201_*(tmp44_);
   double tmp203_ = tmp200_-tmp202_;
   double tmp204_ = 0.500000*(tmp203_);
   double tmp205_ = (tmp204_)/tmp66_;
   double tmp206_ = -(tmp70_);
   double tmp207_ = tmp206_*tmp66_;
   double tmp208_ = (tmp17_)-tmp44_;
   double tmp209_ = -(tmp72_);
   double tmp210_ = tmp209_*tmp66_;
   double tmp211_ = (tmp20_)-tmp44_;
   double tmp212_ = (tmp23_)-tmp44_;
   double tmp213_ = (tmp26_)-tmp44_;
   double tmp214_ = (tmp29_)-tmp44_;
   double tmp215_ = (tmp32_)-tmp44_;
   double tmp216_ = (tmp35_)-tmp44_;

  mVal[0] = (mLocGr1_0-tmp8_)/tmp109_-(tmp67_)/tmp66_;

  mCompDer[0][0] = -(((mLocDGr2Dz0-tmp70_)*tmp66_-(tmp67_)*(tmp119_))/tmp73_);
  mCompDer[0][1] = -((tmp210_-(tmp67_)*(tmp128_))/tmp73_);
  mCompDer[0][2] = -((tmp130_-(tmp67_)*(tmp139_))/tmp73_);
  mCompDer[0][3] = -((tmp141_-(tmp67_)*(tmp150_))/tmp73_);
  mCompDer[0][4] = -((tmp152_-(tmp67_)*(tmp161_))/tmp73_);
  mCompDer[0][5] = -((tmp163_-(tmp67_)*(tmp172_))/tmp73_);
  mCompDer[0][6] = -((tmp174_-(tmp67_)*(tmp183_))/tmp73_);
  mCompDer[0][7] = -((tmp185_-(tmp67_)*(tmp194_))/tmp73_);
  mCompDer[0][8] = -((tmp196_-(tmp67_)*(tmp205_))/tmp73_);
  mVal[1] = (mLocGr1_1-tmp8_)/tmp109_-(tmp110_)/tmp66_;

  mCompDer[1][0] = -((tmp207_-(tmp110_)*(tmp119_))/tmp73_);
  mCompDer[1][1] = -(((mLocDGr2Dz1-tmp72_)*tmp66_-(tmp110_)*(tmp128_))/tmp73_);
  mCompDer[1][2] = -((tmp130_-(tmp110_)*(tmp139_))/tmp73_);
  mCompDer[1][3] = -((tmp141_-(tmp110_)*(tmp150_))/tmp73_);
  mCompDer[1][4] = -((tmp152_-(tmp110_)*(tmp161_))/tmp73_);
  mCompDer[1][5] = -((tmp163_-(tmp110_)*(tmp172_))/tmp73_);
  mCompDer[1][6] = -((tmp174_-(tmp110_)*(tmp183_))/tmp73_);
  mCompDer[1][7] = -((tmp185_-(tmp110_)*(tmp194_))/tmp73_);
  mCompDer[1][8] = -((tmp196_-(tmp110_)*(tmp205_))/tmp73_);
  mVal[2] = (mLocGr1_2-tmp8_)/tmp109_-(tmp208_)/tmp66_;

  mCompDer[2][0] = -((tmp207_-(tmp208_)*(tmp119_))/tmp73_);
  mCompDer[2][1] = -((tmp210_-(tmp208_)*(tmp128_))/tmp73_);
  mCompDer[2][2] = -(((mLocDGr2Dz2-tmp75_)*tmp66_-(tmp208_)*(tmp139_))/tmp73_);
  mCompDer[2][3] = -((tmp141_-(tmp208_)*(tmp150_))/tmp73_);
  mCompDer[2][4] = -((tmp152_-(tmp208_)*(tmp161_))/tmp73_);
  mCompDer[2][5] = -((tmp163_-(tmp208_)*(tmp172_))/tmp73_);
  mCompDer[2][6] = -((tmp174_-(tmp208_)*(tmp183_))/tmp73_);
  mCompDer[2][7] = -((tmp185_-(tmp208_)*(tmp194_))/tmp73_);
  mCompDer[2][8] = -((tmp196_-(tmp208_)*(tmp205_))/tmp73_);
  mVal[3] = (mLocGr1_3-tmp8_)/tmp109_-(tmp211_)/tmp66_;

  mCompDer[3][0] = -((tmp207_-(tmp211_)*(tmp119_))/tmp73_);
  mCompDer[3][1] = -((tmp210_-(tmp211_)*(tmp128_))/tmp73_);
  mCompDer[3][2] = -((tmp130_-(tmp211_)*(tmp139_))/tmp73_);
  mCompDer[3][3] = -(((mLocDGr2Dz3-tmp77_)*tmp66_-(tmp211_)*(tmp150_))/tmp73_);
  mCompDer[3][4] = -((tmp152_-(tmp211_)*(tmp161_))/tmp73_);
  mCompDer[3][5] = -((tmp163_-(tmp211_)*(tmp172_))/tmp73_);
  mCompDer[3][6] = -((tmp174_-(tmp211_)*(tmp183_))/tmp73_);
  mCompDer[3][7] = -((tmp185_-(tmp211_)*(tmp194_))/tmp73_);
  mCompDer[3][8] = -((tmp196_-(tmp211_)*(tmp205_))/tmp73_);
  mVal[4] = (mLocGr1_4-tmp8_)/tmp109_-(tmp212_)/tmp66_;

  mCompDer[4][0] = -((tmp207_-(tmp212_)*(tmp119_))/tmp73_);
  mCompDer[4][1] = -((tmp210_-(tmp212_)*(tmp128_))/tmp73_);
  mCompDer[4][2] = -((tmp130_-(tmp212_)*(tmp139_))/tmp73_);
  mCompDer[4][3] = -((tmp141_-(tmp212_)*(tmp150_))/tmp73_);
  mCompDer[4][4] = -(((mLocDGr2Dz4-tmp79_)*tmp66_-(tmp212_)*(tmp161_))/tmp73_);
  mCompDer[4][5] = -((tmp163_-(tmp212_)*(tmp172_))/tmp73_);
  mCompDer[4][6] = -((tmp174_-(tmp212_)*(tmp183_))/tmp73_);
  mCompDer[4][7] = -((tmp185_-(tmp212_)*(tmp194_))/tmp73_);
  mCompDer[4][8] = -((tmp196_-(tmp212_)*(tmp205_))/tmp73_);
  mVal[5] = (mLocGr1_5-tmp8_)/tmp109_-(tmp213_)/tmp66_;

  mCompDer[5][0] = -((tmp207_-(tmp213_)*(tmp119_))/tmp73_);
  mCompDer[5][1] = -((tmp210_-(tmp213_)*(tmp128_))/tmp73_);
  mCompDer[5][2] = -((tmp130_-(tmp213_)*(tmp139_))/tmp73_);
  mCompDer[5][3] = -((tmp141_-(tmp213_)*(tmp150_))/tmp73_);
  mCompDer[5][4] = -((tmp152_-(tmp213_)*(tmp161_))/tmp73_);
  mCompDer[5][5] = -(((mLocDGr2Dz5-tmp81_)*tmp66_-(tmp213_)*(tmp172_))/tmp73_);
  mCompDer[5][6] = -((tmp174_-(tmp213_)*(tmp183_))/tmp73_);
  mCompDer[5][7] = -((tmp185_-(tmp213_)*(tmp194_))/tmp73_);
  mCompDer[5][8] = -((tmp196_-(tmp213_)*(tmp205_))/tmp73_);
  mVal[6] = (mLocGr1_6-tmp8_)/tmp109_-(tmp214_)/tmp66_;

  mCompDer[6][0] = -((tmp207_-(tmp214_)*(tmp119_))/tmp73_);
  mCompDer[6][1] = -((tmp210_-(tmp214_)*(tmp128_))/tmp73_);
  mCompDer[6][2] = -((tmp130_-(tmp214_)*(tmp139_))/tmp73_);
  mCompDer[6][3] = -((tmp141_-(tmp214_)*(tmp150_))/tmp73_);
  mCompDer[6][4] = -((tmp152_-(tmp214_)*(tmp161_))/tmp73_);
  mCompDer[6][5] = -((tmp163_-(tmp214_)*(tmp172_))/tmp73_);
  mCompDer[6][6] = -(((mLocDGr2Dz6-tmp83_)*tmp66_-(tmp214_)*(tmp183_))/tmp73_);
  mCompDer[6][7] = -((tmp185_-(tmp214_)*(tmp194_))/tmp73_);
  mCompDer[6][8] = -((tmp196_-(tmp214_)*(tmp205_))/tmp73_);
  mVal[7] = (mLocGr1_7-tmp8_)/tmp109_-(tmp215_)/tmp66_;

  mCompDer[7][0] = -((tmp207_-(tmp215_)*(tmp119_))/tmp73_);
  mCompDer[7][1] = -((tmp210_-(tmp215_)*(tmp128_))/tmp73_);
  mCompDer[7][2] = -((tmp130_-(tmp215_)*(tmp139_))/tmp73_);
  mCompDer[7][3] = -((tmp141_-(tmp215_)*(tmp150_))/tmp73_);
  mCompDer[7][4] = -((tmp152_-(tmp215_)*(tmp161_))/tmp73_);
  mCompDer[7][5] = -((tmp163_-(tmp215_)*(tmp172_))/tmp73_);
  mCompDer[7][6] = -((tmp174_-(tmp215_)*(tmp183_))/tmp73_);
  mCompDer[7][7] = -(((mLocDGr2Dz7-tmp85_)*tmp66_-(tmp215_)*(tmp194_))/tmp73_);
  mCompDer[7][8] = -((tmp196_-(tmp215_)*(tmp205_))/tmp73_);
  mVal[8] = (mLocGr1_8-tmp8_)/tmp109_-(tmp216_)/tmp66_;

  mCompDer[8][0] = -((tmp207_-(tmp216_)*(tmp119_))/tmp73_);
  mCompDer[8][1] = -((tmp210_-(tmp216_)*(tmp128_))/tmp73_);
  mCompDer[8][2] = -((tmp130_-(tmp216_)*(tmp139_))/tmp73_);
  mCompDer[8][3] = -((tmp141_-(tmp216_)*(tmp150_))/tmp73_);
  mCompDer[8][4] = -((tmp152_-(tmp216_)*(tmp161_))/tmp73_);
  mCompDer[8][5] = -((tmp163_-(tmp216_)*(tmp172_))/tmp73_);
  mCompDer[8][6] = -((tmp174_-(tmp216_)*(tmp183_))/tmp73_);
  mCompDer[8][7] = -((tmp185_-(tmp216_)*(tmp194_))/tmp73_);
  mCompDer[8][8] = -(((mLocDGr2Dz8-tmp87_)*tmp66_-(tmp216_)*(tmp205_))/tmp73_);
}


void cEqCorrelGrid_9_Im2Var::ComputeValDerivHessian()
{
  ELISE_ASSERT(false,"Foncteur cEqCorrelGrid_9_Im2Var Has no Der Sec");
}

void cEqCorrelGrid_9_Im2Var::SetDGr2Dz0(double aVal){ mLocDGr2Dz0 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetDGr2Dz1(double aVal){ mLocDGr2Dz1 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetDGr2Dz2(double aVal){ mLocDGr2Dz2 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetDGr2Dz3(double aVal){ mLocDGr2Dz3 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetDGr2Dz4(double aVal){ mLocDGr2Dz4 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetDGr2Dz5(double aVal){ mLocDGr2Dz5 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetDGr2Dz6(double aVal){ mLocDGr2Dz6 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetDGr2Dz7(double aVal){ mLocDGr2Dz7 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetDGr2Dz8(double aVal){ mLocDGr2Dz8 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr1_0(double aVal){ mLocGr1_0 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr1_1(double aVal){ mLocGr1_1 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr1_2(double aVal){ mLocGr1_2 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr1_3(double aVal){ mLocGr1_3 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr1_4(double aVal){ mLocGr1_4 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr1_5(double aVal){ mLocGr1_5 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr1_6(double aVal){ mLocGr1_6 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr1_7(double aVal){ mLocGr1_7 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr1_8(double aVal){ mLocGr1_8 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr2of0_0(double aVal){ mLocGr2of0_0 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr2of0_1(double aVal){ mLocGr2of0_1 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr2of0_2(double aVal){ mLocGr2of0_2 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr2of0_3(double aVal){ mLocGr2of0_3 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr2of0_4(double aVal){ mLocGr2of0_4 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr2of0_5(double aVal){ mLocGr2of0_5 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr2of0_6(double aVal){ mLocGr2of0_6 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr2of0_7(double aVal){ mLocGr2of0_7 = aVal;}
void cEqCorrelGrid_9_Im2Var::SetGr2of0_8(double aVal){ mLocGr2of0_8 = aVal;}



double * cEqCorrelGrid_9_Im2Var::AdrVarLocFromString(const std::string & aName)
{
   if (aName == "DGr2Dz0") return & mLocDGr2Dz0;
   if (aName == "DGr2Dz1") return & mLocDGr2Dz1;
   if (aName == "DGr2Dz2") return & mLocDGr2Dz2;
   if (aName == "DGr2Dz3") return & mLocDGr2Dz3;
   if (aName == "DGr2Dz4") return & mLocDGr2Dz4;
   if (aName == "DGr2Dz5") return & mLocDGr2Dz5;
   if (aName == "DGr2Dz6") return & mLocDGr2Dz6;
   if (aName == "DGr2Dz7") return & mLocDGr2Dz7;
   if (aName == "DGr2Dz8") return & mLocDGr2Dz8;
   if (aName == "Gr1_0") return & mLocGr1_0;
   if (aName == "Gr1_1") return & mLocGr1_1;
   if (aName == "Gr1_2") return & mLocGr1_2;
   if (aName == "Gr1_3") return & mLocGr1_3;
   if (aName == "Gr1_4") return & mLocGr1_4;
   if (aName == "Gr1_5") return & mLocGr1_5;
   if (aName == "Gr1_6") return & mLocGr1_6;
   if (aName == "Gr1_7") return & mLocGr1_7;
   if (aName == "Gr1_8") return & mLocGr1_8;
   if (aName == "Gr2of0_0") return & mLocGr2of0_0;
   if (aName == "Gr2of0_1") return & mLocGr2of0_1;
   if (aName == "Gr2of0_2") return & mLocGr2of0_2;
   if (aName == "Gr2of0_3") return & mLocGr2of0_3;
   if (aName == "Gr2of0_4") return & mLocGr2of0_4;
   if (aName == "Gr2of0_5") return & mLocGr2of0_5;
   if (aName == "Gr2of0_6") return & mLocGr2of0_6;
   if (aName == "Gr2of0_7") return & mLocGr2of0_7;
   if (aName == "Gr2of0_8") return & mLocGr2of0_8;
   return 0;
}


cElCompiledFonc::cAutoAddEntry cEqCorrelGrid_9_Im2Var::mTheAuto("cEqCorrelGrid_9_Im2Var",cEqCorrelGrid_9_Im2Var::Alloc);


cElCompiledFonc *  cEqCorrelGrid_9_Im2Var::Alloc()
{  return new cEqCorrelGrid_9_Im2Var();
}



/*Footer-MicMac-eLiSe-25/06/2007

Ce logiciel est un programme informatique servant à la mise en
correspondances d'images pour la reconstruction du relief.

Ce logiciel est régi par la licence CeCILL-B soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL-B telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL-B, et que vous en avez accepté les
termes.
Footer-MicMac-eLiSe-25/06/2007*/
