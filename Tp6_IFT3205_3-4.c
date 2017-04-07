/*------------------------------------------------------*/
/* Prog    : Tp6_IFT3205.c                              */
/* Auteur  : Francois Poitras et Charles Langlois                                       */
/* Date    : --/--/2010                                 */
/* version :                                            */
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo6.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/
/*------------------------------------------------*/
#define NAME_VISUALISER "./ViewSig.sh "

/*------------------------------------------------*/
/* PROTOTYPE DE FONCTIONS  -----------------------*/
/*------------------------------------------------*/
float hamming(int t) {
    return 0.54 - 0.46 * cos(2*PI *t/512.0);
}
/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/
/*------------------------------------------------*/
int memechose31(int argc, char **argv)
{
  int i,j,k;
  int length;
  char BufSystVisuSig[100];

  //===============================
  //Question 2.1.(a)
  //===============================
  float*  Sign1=LoadSignalDat("Spectre",&length);
  float*  Sign1I=fmatrix_allocate_1d(length);
  float*  Sign1M=fmatrix_allocate_1d(length); 
  FFT1D(Sign1,Sign1I,length);
  ModVct(Sign1M,Sign1,Sign1I,length);
  CenterVct(Sign1M,length);
  SaveSignalDat("Spectre",Sign1M,length);  
   
  //Visu
  strcpy(BufSystVisuSig,NAME_VISUALISER);
  strcat(BufSystVisuSig,"FFT_Spectre.dat&");
  printf(" %s",BufSystVisuSig);
  system(BufSystVisuSig);

  //==End=========================================================
  //==============================================================

  //retour sans probleme
  printf("\n C'est fini ... \n\n");
  return 0;    
}

int memechose33(int argc, char **argv)
{
  int i, j, k = 0;
  int length;
  char BufSystVisuSig[100];

  //===============================
  //Question 2.1.(a)
  //===============================
  float*  Sign1 = LoadSignalDat("Spectre", &length);
  float*  Sign1I = fmatrix_allocate_1d(length);
  float*  Sign1M = fmatrix_allocate_1d(length);

  float *buff_r = fmatrix_allocate_1d(512);
  float *buff_i = fmatrix_allocate_1d(512);
  float *buff_rmoy = fmatrix_allocate_1d(512);
  float *buff_imoy = fmatrix_allocate_1d(512);

  for (i = 0; i < 512; i++)
  {
    buff_r[i] = 0;
    buff_i[i] = 0;
    buff_rmoy[i] = 0;
    buff_imoy[i] = 0;
  }

  for (i = 0; i < length - 256; i += 256)
  {
    ++k;
    for (j = 0; j < 512; j++)
    {
      buff_r[j] = hamming(j) * Sign1[i + j];
      buff_i[j] = 0;
    }
    FFT1D(buff_r, buff_i, 512);
    for (j = 0; j < 512; j++) {
      buff_rmoy[j] += fabs(buff_r[j]);
      buff_imoy[j] += fabs(buff_i[j]);
    }
  }


  printf("\nk=%d", k);
  //moyennation
  for (i = 0; i < 512; i++) {
    buff_rmoy[i] /= (float)k;
    buff_imoy[i] /= (float)k;
  }

  ModVct(buff_r, buff_rmoy, buff_imoy, 512);
  CenterVct(buff_r, 512);
  SaveSignalDat("FFT_SpectremoyHamm",buff_r,512);  

  //Visu
  strcpy(BufSystVisuSig, NAME_VISUALISER);
  strcat(BufSystVisuSig, "FFT_SpectremoyHamm.dat&");
  printf(" %s", BufSystVisuSig);
  system(BufSystVisuSig);

  //==End=========================================================
  //==============================================================

  //retour sans probleme
  printf("\n C'est fini ... \n\n");
  return 0;
}

int memechose32(int argc, char **argv)
{
  int i, j, k = 0;
  int length;
  char BufSystVisuSig[100];

  //===============================
  //Question 2.1.(a)
  //===============================
  float*  Sign1 = LoadSignalDat("Spectre", &length);
  float*  Sign1I = fmatrix_allocate_1d(length);
  float*  Sign1M = fmatrix_allocate_1d(length);

  float *buff_r = fmatrix_allocate_1d(512);
  float *buff_i = fmatrix_allocate_1d(512);
  float *buff_rmoy = fmatrix_allocate_1d(512);
  float *buff_imoy = fmatrix_allocate_1d(512);

  for (i = 0; i < 512; i++)
  {
    buff_r[i] = 0;
    buff_i[i] = 0;
    buff_rmoy[i] = 0;
    buff_imoy[i] = 0;
  }

  for (i = 0; i < length - 256; i += 256)
  {
    ++k;
    for (j = 0; j < 512; j++)
    {
      buff_r[j] = Sign1[i + j];
      buff_i[j] = 0;
    }
    FFT1D(buff_r, buff_i, 512);
    for (j = 0; j < 512; j++) {
      buff_rmoy[j] += fabs(buff_r[j]);
      buff_imoy[j] += fabs(buff_i[j]);
    }
  }


  printf("\nk=%d", k);
  //moyennation
  for (i = 0; i < 512; i++) {
    buff_rmoy[i] /= (float)k;
    buff_imoy[i] /= (float)k;
  }

  ModVct(buff_r, buff_rmoy, buff_imoy, 512);
  CenterVct(buff_r, 512);
  SaveSignalDat("FFT_Spectremoy",buff_r,512);  

  //Visu
  strcpy(BufSystVisuSig, NAME_VISUALISER);
  strcat(BufSystVisuSig, "FFT_Spectremoy.dat&");
  printf(" %s", BufSystVisuSig);
  system(BufSystVisuSig);

  //==End=========================================================
  //==============================================================

  //retour sans probleme
  printf("\n C'est fini ... \n\n");
  return 0;
}


int main(int argc, char **argv)
{
    memechose31(argc, argv);
    memechose32(argc, argv);
    memechose33(argc, argv);
}