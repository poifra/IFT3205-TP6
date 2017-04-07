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

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL   -------------------------*/
/*------------------------------------------------*/
int main(int argc, char **argv)
{
  int i, j, k = 0;
  int length;
  char BufSystVisuSig[100];

  //===============================
  //Question 2.1.(a)
  //===============================
  float*  Sign1 = LoadSignalDat("moteur1", &length);
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
  SaveSignalDat("FFT_Moteur2moy",buff_r,512);  

  //Visu
  strcpy(BufSystVisuSig, NAME_VISUALISER);
  strcat(BufSystVisuSig, "FFT_Moteur2moy.dat&");
  printf(" %s", BufSystVisuSig);
  system(BufSystVisuSig);

  //==End=========================================================
  //==============================================================

  //retour sans probleme
  printf("\n C'est fini ... \n\n");
  return 0;
}


