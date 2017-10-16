/*
Copyright (c) 2015, Patrick Weltevrede
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "psrsalsa.h"







#define PUMA_STD 






#ifndef _PUMA_H
#define _PUMA_H 
#define MMAX(A,B) (A > B? A: B)
#define MAXDMS 4
#define MAXSTOKESPARAMS 4
#define MAXPACKS (MMAX(MAXDMS, MAXSTOKESPARAMS))
#define PACKX 0
#define PACKY 1
#define PACKI 0
#define PACKQ 1
#define PACKU 2
#define PACKV 3
#define MAXFFTSIZE 8192
#define MINFFTSIZE 32
#define NoMODE -9
#define MODEM1 -1
#define MODE0 0
#define MODE1 1
#define MODE2 2
#define MODE3 3
#define MODE4 4
#ifndef TIME
#include <time.h>
#ifndef Boolean
#define Boolean int
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#define NOREQUESTRECEIVED 0
#define WORKINGONREQUEST 1
#define IGNORE "_IgNoRe_"
#define TXTLEN 64
#define NAMELEN 16
#define COMMENTLEN 256
#define FILELEN 64
#define CMDLEN 256
#define DIRLEN 256
#define PATHLEN 256
#define NCRATES 2
#define CRATE0 0
#define CRATE1 1
#define MODE0REP "mode0rep"
#define MODE1REP "mode1rep"
#define ADJUST "adjust"
#define ARCREP "arcrep"
#define DELREP "delrep"
#define BINDAT "bindat"
#define OBSFNS "obsfns"
#define ASTRO "puma"
#define NOFCLUSTERS 4
#define MAXDISKS 4
#define MAXPARTITIONS 5
#define MAXDSPSPERBOARD 6
#define MAXBOARDSPERCLUSTER 4
#define MAXDSPBOARDS 16
#define MAXDSPS (NOFCLUSTERS * MAXBOARDSPERCLUSTER * MAXDSPSPERBOARD)
#define MAXFREQBANDS 8
#define MAXTELESCOPES 14
#define NoRequest 0
#define GoOffLine -1
#define Initialise 1
#define StartObservation 2
#define AbortObservation 3
#define GetNofFrames 4
#define GetDiskInfo 5
#define GetCrateNumber 6
#define CopyToArchive 7
#define AbortCopy 8
#define DeleteFiles 9
#define GetGeneralInfo 1
#define GetTime 2
#define WantOffLine 3
#define NewTape 4
#define NoError 0x00000000
#define AnyError 0x00000001
#define DSPDown 0x00000002
#define DSPBoardDown 0x00000004
#define ClusterDown 0x00000008
#define AnalogChannelDown 0x00000010
#define TemperatureWarning 0x00000020
#define PowerSupplyWarning 0x00000040
#define LastRequestFailed 0x00000080
#define DiskFull 0x00000100
#define OutOfSync 0x00000200
#define ClashesHappened 0x00000400
#define HPMissedFrame 0x00000800
#define WrongScanObsNumber 0x00001000
#define DeleteFailed 0x00004000
#define DataCorrupt 0x00008000
#define SHARC_Message 0x00010000
#define KilledBySignal 0x00400000
#define TapeError 0x00800000
#define LongTapeDelay 0x01000000
#define CopyAborted 0x02000000
#define ObsAborted 0x04000000
#define ForcedToAbort 0x08000000
#define RevisionError 0x10000000
#define SoftwareError 0x20000000
#define ExceptionReport 0x40000000
#define UnknownError 0x80000000
#define PuMaBootAndTest -1
#define PuMaIdle 0
#define PuMaInitialising 1
#define PuMaAwaitStart 2
#define PuMaObserving 4
#define PuMaAborting 8
#define PuMaArchiving 16
#define PuMaDeleting 32
typedef int RequestType;
typedef int RequestStatus;
typedef long int MJD_t;
typedef char ProgNameType[PATHLEN];
typedef enum {Observation, Test, Calibration} ObsReason;
typedef struct { int i0, i1;} ExtInt;
typedef struct { double itsDeclination;
   double itsRightAscension;
               } SkyPosType;
typedef struct { char itsName[NAMELEN];
                 ObsReason itsReason;
                 SkyPosType itsSkyPosition;
                 char itsEpoch[NAMELEN];
               } SourceType;
typedef struct { char itsObsID [NAMELEN];
   char itsCalObsID[NAMELEN];
                 int itsModeNo;
   MJD_t itsStartDate;
   PUMA_STD time_t itsStartTime;
   long int itsDuration;
   SourceType itsSource;
   char itsPI[TXTLEN];
   char itsProposalNumber[TXTLEN];
        } ObservationType;
typedef struct { Boolean isGeografical;
   double itsLongitude;
   double itsLatitude;
        } CoordType;
typedef struct { char itsName[TXTLEN];
   CoordType itsCoord;
   double itsHeight;
   MJD_t itsCurrentDate;
   PUMA_STD time_t itsCurrentTime;
   int itsGPS_MASERoffset;
   Boolean isUniverseClosed;
               } ObservatoryType;
typedef struct { Boolean isAvailable;
                 char itsFrontendID[TXTLEN];
   char itsFrontendStatus[TXTLEN];
   Boolean isNoiseSrcOff;
   double itsSysTemp;
               } TelescopeType;
typedef struct { double itsMidFreq;
                 double itsWidth;
                 double itsMidSkyFreq;
                 Boolean isNeeded;
                 Boolean isFlipped;
        } BandType;
typedef struct { char itsDir[DIRLEN];
   char itsFile[FILELEN];
   char itsRevision[TXTLEN];
        } SWType;
typedef struct { Boolean isOffsetDynamic;
   Boolean isScaleDynamic;
   long int itsAdjustInterval;
        } AdjustType;
typedef struct { char itsComment[COMMENTLEN];
   char itsCrateNumber;
   ObservatoryType itsObservatory;
   ObservationType itsReqObs;
   TelescopeType itsTelescopes[MAXTELESCOPES];
   BandType itsFreqBands [MAXFREQBANDS];
   int itsBandToClusMap[MAXFREQBANDS];
   AdjustType itsAdjust;
   SWType itsDSPTestProg;
   char itsWarnDigiTemp;
   int itsDATTapeSize;
               } InitType;
typedef struct { Boolean isXUsed;
                 Boolean isYUsed;
               } PolUsedType;
typedef struct { SWType itsSW;
   int itsDecimation;
        } FilterType;
typedef struct { InitType itsInit;
                 SWType itsMasterSoftware;
   SWType itsDSPSoftware [NOFCLUSTERS][MAXDSPSPERBOARD];
   float itsScaleFactor [NOFCLUSTERS][MAXPACKS];
   float itsDCvalue [NOFCLUSTERS][MAXPACKS];
   float itsXPolScaleFac[NOFCLUSTERS];
   double itsNofSigmas;
                 PolUsedType itsBandPol;
                 int itsNofBitsPerSample;
                 FilterType itsFilter;
                 Boolean isStokesParamNeeded[MAXSTOKESPARAMS];
                 int itsNofFFTFreqChannels;
                 int itsNofSampsToAdd;
                 int itsNofSHARCsToAdd;
                 double itsDM[MAXDMS];
   double itsRM;
               } ModeType;
typedef struct { char itsPath[PATHLEN];
               } TaskType;
typedef struct { long int itsMaxSpace [MAXPARTITIONS];
                 long int itsFreeSpace[MAXPARTITIONS];
               } DiskType;
typedef struct { Boolean isPos5VoltOK;
   Boolean isNeg5VoltOK;
        } AnaPwrType;
typedef struct { Boolean isPos5VoltOK;
        } DigPwrType;
typedef struct { Boolean isClusterOK;
                 Boolean isAnalogOK;
                 Boolean isAnalogTempOK;
                 AnaPwrType itsAnaPwrSupply;
               } ClusterWarnType;
typedef struct { long int itsNofClashes;
                 Boolean isDSPBoardOK;
                 Boolean isDSPBoardTempOK;
                 DigPwrType itsDigPwrSupply;
               } DSPBoardWarnType;
typedef struct { long int itsNofObsFiles;
   char itsDir[DIRLEN];
               } StoreLocType;
typedef struct { char itsCrateNumber;
   volatile long int itsState;
   unsigned long int itsErrorMessage;
   long int itsNofCounter;
   long int itsTotCounter;
   char itsUnit[TXTLEN];
   ClusterWarnType itsCluster [NOFCLUSTERS];
   char itsADCNumber[NOFCLUSTERS];
   DSPBoardWarnType itsDSPBoard [MAXDSPBOARDS];
   ExtInt itsBoardNumber[MAXDSPBOARDS];
   DiskType itsDisks [MAXDISKS];
   PUMA_STD time_t itsActualStartTime;
   StoreLocType itsLastObsFiles;
   Boolean isDSPOK[MAXDSPS];
   float itsLastScaleFactor[NOFCLUSTERS][MAXPACKS];
   float itsLastDCValue [NOFCLUSTERS][MAXPACKS];
   float itsLastXYRatio [NOFCLUSTERS];
   int itsLastModeNumber;
        } StatusType;
typedef struct { volatile RequestType itsRequestTMStoPuMa;
     RequestStatus itsRequestState;
     ModeType itsMode;
     TaskType itsOtherTask;
     ProgNameType itsProgName;
                 volatile Boolean isActive;
               } TMSAreaType;
typedef struct { volatile RequestType itsRequestPuMatoTMS;
                          RequestStatus itsRequestState;
                          StatusType itsPuMaStatus;
                          ProgNameType itsProgName;
   volatile Boolean isActive;
               } PuMaAreaType;
#define TMSSHMNAME "/shmemtms"
#define PuMaSHMNAME "/shmempuma"
#endif
#endif
#ifndef _PUMADATA_H
#define _PUMADATA_H 
#define HEADERFORMATID "DPC_1.3"
#define defRednSoftwareVer 1
#define ASTROFILESIZE (200 * 1024 * 1024)
#define LONGNAMELEN 24
typedef char Dummy;
typedef struct { Boolean Active;
     Boolean NoiseSrcOff;
     char FrontendID[TXTLEN];
     char FrontendStatus[TXTLEN];
     double Tsys;
        } Telescope_type;
typedef struct { Boolean Operational;
                   Boolean NonFlip;
                   double MidFreq;
                   double Width;
                   double SkyMidFreq;
               } Band_type;
typedef struct { char HdrVer[NAMELEN];
     char Platform[NAMELEN];
     char ThisFileName[LONGNAMELEN];
     char ScanNum[NAMELEN];
     char Comment[COMMENTLEN];
     int NFiles;
     int FileNum;
     char TapeID[LONGNAMELEN];
     int NTapes;
     int TapeNum;
     int ParBlkSize;
     int DataBlkSize;
     Boolean Cluster[MAXFREQBANDS];
     int DataMJD;
                   Dummy AlignDouble[4];
     double DataTime;
     Dummy GeneralDummy[64];
        } General_type;
typedef struct { char Name[NAMELEN];
                   double Long;
                   double Lat;
                   double Height;
               } Observatory_type;
typedef struct { char ObsName[TXTLEN];
                   char ObsType[TXTLEN];
                   char LastCalID[NAMELEN];
                   int StMJD;
                   int StTime;
                   double StLST;
                   double MaserOffset;
     double Dur;
                   double IonRM;
                   Dummy ObservationDummy[64];
               } Observation_type;
typedef struct { char Pulsar[NAMELEN];
                   double RA;
                   double Dec;
                   char Epoch[NAMELEN];
               } Target_type;
typedef struct { Telescope_type Tel[MAXTELESCOPES];
     Band_type Band[MAXFREQBANDS];
                   char Backend[NAMELEN];
                   Boolean ResFringeOff;
                   Boolean AddBoxOn;
     int BandsToClusMap[MAXFREQBANDS];
                   Dummy WSRTDummy[16];
        } Signalpath_type;
typedef struct { int Nr;
                   Dummy SpareDummy1[4];
     float XPolScaleFac[MAXFREQBANDS];
                   Boolean ActCluster[MAXFREQBANDS];
     int FIRFactor;
     int NSHARCsAdded;
     int NSampsAdded;
     int NFreqInFile;
                   int Tsamp;
                   Boolean Iout;
                   Boolean Qout;
                   Boolean Uout;
                   Boolean Vout;
                   Boolean Xout;
                   Boolean Yout;
                   Dummy SpareDummy2[8];
                   int NDMs;
                   Dummy SpareDummy3[4];
                   Dummy AlignDouble[4];
     double DM[MAXDMS];
                   double RM;
                   Boolean DC_Dynamic;
                   Boolean ScaleDynamic;
                   double AdjustInterval;
                   Boolean FloatsOut;
                   int BitsPerSamp;
     double SigmaRange;
                   Dummy ModeDummy[64 - 8];
        } Mode_type;
typedef struct { char Master[TXTLEN];
                   char DSP[MAXDSPSPERBOARD][TXTLEN];
     char Filter[TXTLEN];
        } Software_type;
typedef struct { int ExitCodeObsy;
     int ExitCodePuMa;
     int ExitCodeClust[MAXFREQBANDS];
     int ExitCodeDataConsistency;
     Dummy CheckDummy[20];
        } Check_type;
typedef struct { Boolean Raw;
     int MJDint;
     double MJDfrac;
     double DM;
     int NTimeInts;
     int NFreqs;
     double DeltaTime;
     double FreqCent;
     double DeltaFreq;
     Boolean IsDedisp;
                   Boolean IsCohDedisp;
                   int CohFFTSize;
     Boolean IsPwr;
     char Zapfile[NAMELEN];
     int NBins;
     Boolean Folded;
     double FoldPeriod;
     Boolean Polyco;
     int NCoef;
     int PolycoSpan;
                   Boolean IsAdjusted;
                   Boolean PolycoStored;
     Boolean Bary;
     Boolean OI;
     Boolean OQ;
     Boolean OU;
     Boolean OV;
     Boolean OP;
     Boolean OTheta;
     Boolean Op;
     Boolean Ov;
     Boolean Opoldeg;
     Boolean OX;
     Boolean OY;
     int TRedn;
                   char Command[CMDLEN];
     int RednSoftwareVer;
                   int GenType;
                  Dummy ReductionDummy[56];
        } Reduction_type;
#define PuMa_Standard_Profile 0x0001
#define PuMa_Autocorrelation_Spectrum 0x0002
#define PuMa_Gated_Data 0x0004
#define PuMa_Dynamic_Spectrum 0x0008
#define PuMa_Reference_Flux_Density 0x0010
#define PuMa_Absolute_Flux_Density 0x0020
#define PuMa_Feed_Corrected 0x0040
#define PuMa_Platform_Corrected 0x0080
#define PuMa_Polarization_Calibrated 0x0100
#define PuMa_RM_Corrected 0x0200
typedef struct { General_type gen;
     Observatory_type obsy;
     Observation_type obs;
     Target_type src;
     Signalpath_type WSRT;
     Mode_type mode;
     Software_type software;
     Check_type check;
     Reduction_type redn;
        } Header_type;
typedef struct { int framenumber;
   int whatfor;
   float scale;
   float offset;
        } Adjustments;
#endif
int GetBEint(int *src)
{
  int dst;
  ((char*)&dst)[0] = ((char*)src)[3];
  ((char*)&dst)[1] = ((char*)src)[2];
  ((char*)&dst)[2] = ((char*)src)[1];
  ((char*)&dst)[3] = ((char*)src)[0];
  return dst;
}
short GetBEshort(short *src)
{
  short dst;
  ((char*)&dst)[0] = ((char*)src)[1];
  ((char*)&dst)[1] = ((char*)src)[0];
  return dst;
}
float GetBEfloat(float *src)
{
  float dst;
  ((char*)&dst)[0] = ((char*)src)[3];
  ((char*)&dst)[1] = ((char*)src)[2];
  ((char*)&dst)[2] = ((char*)src)[1];
  ((char*)&dst)[3] = ((char*)src)[0];
  return dst;
}
double GetBEdouble(double *src)
{
  double dst;
  ((char*)&dst)[0] = ((char*)src)[7];
  ((char*)&dst)[1] = ((char*)src)[6];
  ((char*)&dst)[2] = ((char*)src)[5];
  ((char*)&dst)[3] = ((char*)src)[4];
  ((char*)&dst)[4] = ((char*)src)[3];
  ((char*)&dst)[5] = ((char*)src)[2];
  ((char*)&dst)[6] = ((char*)src)[1];
  ((char*)&dst)[7] = ((char*)src)[0];
  return dst;
}
void PutBEint(int *dst, int src)
{
  ((char*)dst)[0] = ((char*)&src)[3];
  ((char*)dst)[1] = ((char*)&src)[2];
  ((char*)dst)[2] = ((char*)&src)[1];
  ((char*)dst)[3] = ((char*)&src)[0];
}
void PutBEshort(short *dst, short src)
{
  ((char*)dst)[0] = ((char*)&src)[1];
  ((char*)dst)[1] = ((char*)&src)[0];
}
void PutBEfloat(float *dst, float src)
{
  ((char*)dst)[0] = ((char*)&src)[3];
  ((char*)dst)[1] = ((char*)&src)[2];
  ((char*)dst)[2] = ((char*)&src)[1];
  ((char*)dst)[3] = ((char*)&src)[0];
}
void PutBEdouble(double *dst, double src)
{
  ((char*)dst)[0] = ((char*)&src)[7];
  ((char*)dst)[1] = ((char*)&src)[6];
  ((char*)dst)[2] = ((char*)&src)[5];
  ((char*)dst)[3] = ((char*)&src)[4];
  ((char*)dst)[4] = ((char*)&src)[3];
  ((char*)dst)[5] = ((char*)&src)[2];
  ((char*)dst)[6] = ((char*)&src)[1];
  ((char*)dst)[7] = ((char*)&src)[0];
}
void swapWord(void *xptr, int size)
{
  unsigned char *cptr = (unsigned char *)xptr;
  unsigned char tmp;
  int i;
  for (i=0; i < size/2; i++)
  {
    tmp = cptr[i];
    cptr[i] = cptr[size-1-i];
    cptr[size-1-i] = tmp;
  }
}
void swapArray(void *xptr, int nelem, int size)
{
  int i;
  unsigned char *cptr = (unsigned char *)xptr;
  for (i=0; i < nelem; i++)
  {
    swapWord(cptr, size);
    cptr += size;
  }
}
void beheader_convert(Header_type *src, Header_type *dst)
{
  int i;
  *dst = *src;
    for (i=0; i< MAXTELESCOPES; i++) {
      dst->WSRT.Tel[i].Active = GetBEint(&src->WSRT.Tel[i].Active);
      dst->WSRT.Tel[i].NoiseSrcOff = GetBEint(&src->WSRT.Tel[i].NoiseSrcOff);
      dst->WSRT.Tel[i].Tsys = GetBEdouble(&src->WSRT.Tel[i].Tsys);
  }
    for (i=0; i< MAXFREQBANDS; i++) {
      dst->WSRT.Band[i].Operational = GetBEint(&src->WSRT.Band[i].Operational);
      dst->WSRT.Band[i].NonFlip = GetBEint(&src->WSRT.Band[i].NonFlip);
      dst->WSRT.BandsToClusMap[i] = GetBEint(&src->WSRT.BandsToClusMap[i]);
      dst->gen.Cluster[i] = GetBEint(&src->gen.Cluster[i]);
      dst->WSRT.Band[i].MidFreq = GetBEdouble(&src->WSRT.Band[i].MidFreq);
      dst->WSRT.Band[i].Width = GetBEdouble(&src->WSRT.Band[i].Width);
      dst->WSRT.Band[i].SkyMidFreq = GetBEdouble(&src->WSRT.Band[i].SkyMidFreq);
    }
    dst->gen.NFiles = GetBEint(&src->gen.NFiles);
    dst->gen.FileNum = GetBEint(&src->gen.FileNum);
    dst->gen.NTapes = GetBEint(&src->gen.NTapes);
    dst->gen.TapeNum = GetBEint(&src->gen.TapeNum);
    dst->gen.ParBlkSize = GetBEint(&src->gen.ParBlkSize);
    dst->gen.DataBlkSize = GetBEint(&src->gen.DataBlkSize);
    dst->gen.DataMJD = GetBEint(&src->gen.DataMJD);
    dst->gen.DataTime = GetBEdouble(&src->gen.DataTime);
    dst->obsy.Long = GetBEdouble(&src->obsy.Long);
    dst->obsy.Lat = GetBEdouble(&src->obsy.Lat);
    dst->obsy.Height = GetBEdouble(&src->obsy.Height);
    dst->obs.StMJD = GetBEint(&src->obs.StMJD);
    dst->obs.StTime = GetBEint(&src->obs.StTime);
    dst->obs.StLST = GetBEdouble(&src->obs.StLST);
    dst->obs.MaserOffset = GetBEdouble(&src->obs.MaserOffset);
    dst->obs.Dur = GetBEdouble(&src->obs.Dur);
    dst->obs.IonRM = GetBEdouble(&src->obs.IonRM);
    dst->src.RA = GetBEdouble(&src->src.RA);
    dst->src.Dec = GetBEdouble(&src->src.Dec);
    dst->WSRT.ResFringeOff = GetBEint(&src->WSRT.ResFringeOff);
    dst->WSRT.AddBoxOn = GetBEint(&src->WSRT.AddBoxOn);
    dst->mode.Nr = GetBEint(&src->mode.Nr);
    dst->mode.FIRFactor = GetBEint(&src->mode.FIRFactor);
    dst->mode. NSHARCsAdded= GetBEint(&src->mode.NSHARCsAdded);
    dst->mode.NSampsAdded = GetBEint(&src->mode.NSampsAdded);
    dst->mode.NFreqInFile = GetBEint(&src->mode.NFreqInFile);
    dst->mode.Tsamp = GetBEint(&src->mode.Tsamp);
    dst->mode.Iout = GetBEint(&src->mode.Iout);
    dst->mode.Qout = GetBEint(&src->mode.Qout);
    dst->mode.Uout = GetBEint(&src->mode.Uout);
    dst->mode.Vout = GetBEint(&src->mode.Vout);
    dst->mode.Xout = GetBEint(&src->mode.Xout);
    dst->mode.Yout = GetBEint(&src->mode.Yout);
    dst->mode.NDMs = GetBEint(&src->mode.NDMs);
    dst->mode.DC_Dynamic = GetBEint(&src->mode.DC_Dynamic);
    dst->mode.ScaleDynamic = GetBEint(&src->mode.ScaleDynamic);
    dst->mode.FloatsOut = GetBEint(&src->mode.FloatsOut);
    dst->mode.BitsPerSamp = GetBEint(&src->mode.BitsPerSamp);
    for (i=0;i<MAXDMS;i++)
      dst->mode.DM[i] = GetBEdouble(&src->mode.DM[i]);
    dst->mode.RM = GetBEdouble(&src->mode.RM);
    dst->mode.AdjustInterval = GetBEdouble(&src->mode.AdjustInterval);
    dst->mode.SigmaRange = GetBEdouble(&src->mode.SigmaRange);
    for (i=0;i<MAXFREQBANDS;i++)
      {
 dst->mode.XPolScaleFac[i] = GetBEfloat(&src->mode.XPolScaleFac[i]);
 dst->mode.ActCluster[i] = GetBEint(&src->mode.ActCluster[i]);
 dst->check.ExitCodeClust[i] = GetBEint(&src->check.ExitCodeClust[i]);
      }
    dst->check.ExitCodeObsy = GetBEint(&src->check.ExitCodeObsy);
    dst->check.ExitCodePuMa = GetBEint(&src->check.ExitCodePuMa);
    dst->check.ExitCodeDataConsistency = GetBEint(&src->check.ExitCodeDataConsistency);
    dst->redn.Raw = GetBEint(&src->redn.Raw);
    dst->redn.MJDint = GetBEint(&src->redn.MJDint);
    dst->redn.NTimeInts = GetBEint(&src->redn.NTimeInts);
    dst->redn.NFreqs = GetBEint(&src->redn.NFreqs);
    dst->redn.IsDedisp = GetBEint(&src->redn.IsDedisp);
    dst->redn.IsCohDedisp = GetBEint(&src->redn.IsCohDedisp);
    dst->redn.CohFFTSize = GetBEint(&src->redn.CohFFTSize);
    dst->redn.IsPwr = GetBEint(&src->redn.IsPwr);
    dst->redn.NBins = GetBEint(&src->redn.NBins);
    dst->redn.Folded = GetBEint(&src->redn.Folded);
    dst->redn.Polyco = GetBEint(&src->redn.Polyco);
    dst->redn.NCoef = GetBEint(&src->redn.NCoef);
    dst->redn.PolycoSpan = GetBEint(&src->redn.PolycoSpan);
    dst->redn.IsAdjusted = GetBEint(&src->redn.IsAdjusted);
    dst->redn.PolycoStored = GetBEint(&src->redn.PolycoStored);
    dst->redn.Bary = GetBEint(&src->redn.Bary);
    dst->redn.OI = GetBEint(&src->redn.OI);
    dst->redn.OQ = GetBEint(&src->redn.OQ);
    dst->redn.OU = GetBEint(&src->redn.OU);
    dst->redn.OV = GetBEint(&src->redn.OV);
    dst->redn.OP = GetBEint(&src->redn.OP);
    dst->redn.OTheta = GetBEint(&src->redn.OTheta);
    dst->redn.Op = GetBEint(&src->redn.Op);
    dst->redn.Ov = GetBEint(&src->redn.Ov);
    dst->redn.Opoldeg = GetBEint(&src->redn.Opoldeg);
    dst->redn.OX = GetBEint(&src->redn.OX);
    dst->redn.OY = GetBEint(&src->redn.OY);
    dst->redn.TRedn = GetBEint(&src->redn.TRedn);
    dst->redn.MJDfrac = GetBEdouble(&src->redn.MJDfrac);
    dst->redn.DM = GetBEdouble(&src->redn.DM);
    dst->redn.DeltaTime = GetBEdouble(&src->redn.DeltaTime);
    dst->redn.FreqCent = GetBEdouble(&src->redn.FreqCent);
    dst->redn.DeltaFreq = GetBEdouble(&src->redn.DeltaFreq);
    dst->redn.FoldPeriod = GetBEdouble(&src->redn.FoldPeriod);
}
int prheader(Header_type *inphdr,FILE *srcfile)
{
  int r;
  Header_type behdr;
#if defined ( __linux__) || defined (__alpha)
  r = fread(&behdr,sizeof(Header_type),1,srcfile);
  beheader_convert(&behdr,inphdr);
#else
  r = fread(inphdr,sizeof(Header_type),1,srcfile);
#endif
  if (strcmp(inphdr->gen.HdrVer,"DPC_1.1")==0 || strcmp(inphdr->gen.HdrVer,"DPC_1.2")==0)
    {
         inphdr->redn.RednSoftwareVer = 0;
    }
   else
   {
       inphdr->redn.RednSoftwareVer = defRednSoftwareVer;
   }
  return r;
}
int pwheader(Header_type *outphdr,FILE *dstfile)
{
  int r;
  Header_type behdr;
#if defined ( __linux__) || defined (__alpha)
  beheader_convert(outphdr,&behdr);
  r = fwrite(&behdr,sizeof(Header_type),1,dstfile);
# else
  r = fwrite(outphdr,sizeof(Header_type),1,dstfile);
#endif
  return r;
}
int puma_nrpol(Header_type hdr)
{
  int NrPol;
  NrPol = 0;
  if((hdr.redn.OX) && (!hdr.redn.OY)) NrPol=1;
  if((!hdr.redn.OX) && (hdr.redn.OY)) NrPol=1;
  if((hdr.redn.OX) && (hdr.redn.OY)) NrPol=2;
  if(hdr.redn.OI) {
      NrPol = 1;
      if(hdr.redn.OQ && hdr.redn.OU && hdr.redn.OV) NrPol=4;
      if(hdr.redn.OP && hdr.redn.OV && hdr.redn.OTheta && hdr.redn.Op && hdr.redn.Ov == 0 && hdr.redn.Opoldeg == 0 && hdr.redn.OU == 0 && hdr.redn.OX == 0 && hdr.redn.OY == 0) NrPol=5;
      if(hdr.redn.OP && hdr.redn.OV && hdr.redn.OTheta && hdr.redn.Op && hdr.redn.Ov == 0 && hdr.redn.Opoldeg == 0 && hdr.redn.OU == 1 && hdr.redn.OX == 1 && hdr.redn.OY == 1) NrPol=8;
  }
  return NrPol;
}
#define ConvertArrayFromBE(xptr,n,size) swapArray(xptr, n, size)
void pumaread(void *prptr, int size, int nelem,FILE *in)
{
#ifdef __alpha
  fread(prptr,size,nelem,in);
  ConvertArrayFromBE(prptr,nelem,size);
#endif
#ifdef __linux__
  fread(prptr,size,nelem,in);
  ConvertArrayFromBE(prptr,nelem,size);
#endif
#ifdef __hpux
  fread(prptr,size,nelem,in);
#endif
}
int pumawrite(void *prptr, int size, int nelem,FILE *out)
{
  int status;
#ifdef __linux__
  ConvertArrayFromBE(prptr,nelem,size);
  status = fwrite(prptr,size,nelem,out);
  ConvertArrayFromBE(prptr,nelem,size);
# else
  status = fwrite(prptr,size,nelem,out);
# endif
  return status;
}
int readWSRTHeader(datafile_definition *datafile, verbose_definition verbose)
{
  int i, version;
  char *dummy, bytevalue;
  int *idummy;
  float *fdummy;
  Header_type puma_hdr;
  prheader(&puma_hdr, datafile->fptr_hdr);
  datafile->NrSubints = puma_hdr.redn.NTimeInts;
  datafile->NrBins = puma_hdr.redn.NBins;
  datafile->NrBits = 8*sizeof(float);
  if(puma_hdr.redn.Folded == 0) {
    if(verbose.debug) {
      printf("DEBUG: redn.Folded flag is not set, assume data is searchmode data\n");
    }
    datafile->NrSubints = 1;
    datafile->gentype = GENTYPE_SEARCHMODE;
    datafile->isFolded = 0;
    datafile->foldMode = FOLDMODE_UNKNOWN;
    datafile->fixedPeriod = -1;
  }else {
    if(verbose.debug) {
      printf("DEBUG: redn.Folded flag is set\n");
    }
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    datafile->fixedPeriod = puma_hdr.redn.FoldPeriod;
  }
  datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  datafile->fixedtsamp = puma_hdr.redn.DeltaTime;
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Memory allocation error");
    return 0;
  }
  datafile->tsub_list[0] = puma_hdr.obs.Dur/(double)datafile->NrSubints;
  datafile->mjd_start = puma_hdr.redn.MJDint + puma_hdr.redn.MJDfrac;
  datafile->freqMode = FREQMODE_UNIFORM;
  if(datafile->freqlabel_list != NULL) {
    free(datafile->freqlabel_list);
    datafile->freqlabel_list = NULL;
  }
  set_centre_frequency(datafile, puma_hdr.redn.FreqCent, verbose);
  double bw;
  bw = puma_hdr.WSRT.Band[0].Width;
  double channelbw = puma_hdr.redn.DeltaFreq;
  if(channelbw < 0 || bw < 0) {
    bw = -fabs(bw);
  }
  if(set_bandwidth(datafile, bw, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Bandwidth changing failed.");
    return 0;
  }
  datafile->NrFreqChan = puma_hdr.redn.NFreqs;
  datafile->ra = puma_hdr.src.RA;
  if(datafile->ra < 0)
    datafile->ra += 2*M_PI;
  datafile->dec = puma_hdr.src.Dec;
  datafile->dm = puma_hdr.redn.DM;
  datafile->rm = puma_hdr.mode.RM;
  datafile->isDeDisp = puma_hdr.redn.IsDedisp;
  if(set_psrname_PSRData(datafile, puma_hdr.src.Pulsar, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting pulsar name failed.");
    return 0;
  }
  if(set_observatory_PSRData(datafile, puma_hdr.obsy.Name, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting observatory name failed.");
    return 0;
  }
  if(verbose.debug) {
    printf("DEBUG: defined long=%lf deg lat=%lf deg height=%lf m\n", puma_hdr.obsy.Long*180.0/M_PI, puma_hdr.obsy.Lat*180.0/M_PI, puma_hdr.obsy.Height);
  }
  tempo2_GRS80_to_ITRF(puma_hdr.obsy.Long, puma_hdr.obsy.Lat, puma_hdr.obsy.Height, &(datafile->telescope_X), &(datafile->telescope_Y), &(datafile->telescope_Z));
  if(verbose.debug) {
    printf("DEBUG: derived ITRF (X,Y,Z) m = (%lf, %lf, %lf) m\n", datafile->telescope_X, datafile->telescope_Y, datafile->telescope_Z);
  }
  double telc_lat, telc_long;
  telc_long = observatory_long_geodetic(*datafile);
  telc_lat = observatory_lat_geodetic(*datafile);
  if(verbose.debug) {
    printf("DEBUG: long=%lf deg lat=%lf (should be the header defined position)\n", telc_long*180.0/M_PI, telc_lat*180.0/M_PI);
  }
  if(telc_long*180.0/M_PI < 6.61 && telc_long*180.0/M_PI > 6.60) {
    if(telc_lat*180.0/M_PI < 53.0 && telc_lat*180.0/M_PI > 52.7) {
      datafile->telescope_X = 3828445.659;
      datafile->telescope_Y = 445223.600000;
      datafile->telescope_Z = 5064921.5677;
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readWSRTHeader: Updated derived ITRF position of telescope with accurate WSRT position.");
    }
  }
  if(set_institute_PSRData(datafile, puma_hdr.gen.Comment, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting institute name failed.");
    return 0;
  }
  for(i = 0; i < strlen(datafile->psrname); i++) {
    if(datafile->psrname[i] == ' ')
      datafile->psrname[i] = '_';
  }
  for(i = 0; i < strlen(datafile->institute); i++) {
    if(datafile->institute[i] == ' ')
      datafile->institute[i] = '_';
  }
  if(set_institute_PSRData(datafile, "Converted", verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting institute name failed.");
    return 0;
  }
  if(set_instrument_PSRData(datafile, "PuMa", verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting instrument name failed.");
    return 0;
  }
  datafile->NrPols = puma_nrpol(puma_hdr);
  if(puma_hdr.redn.OI == TRUE && puma_hdr.redn.OP == TRUE && puma_hdr.redn.OV == TRUE && puma_hdr.redn.OQ == FALSE && puma_hdr.redn.OU == FALSE && puma_hdr.redn.OTheta == FALSE && puma_hdr.redn.Op == FALSE && puma_hdr.redn.Ov == FALSE && puma_hdr.redn.Opoldeg == FALSE && puma_hdr.redn.OX == FALSE && puma_hdr.redn.OY == FALSE) {
    datafile->NrPols = 3;
  }
  dummy = puma_hdr.redn.ReductionDummy;
  i = 1;
  if(dummy[0] != 'p')
    i = 0;
  dummy ++;
  if(dummy[0] != 's')
    i = 0;
  dummy ++;
  if(dummy[0] != 'o')
    i = 0;
  dummy ++;
  if(dummy[0] != 'f')
    i = 0;
  dummy ++;
  if(dummy[0] != 't')
    i = 0;
  dummy ++;
  if(i == 1) {
    idummy = (int *)dummy;
    dummy += sizeof(int);
    version = *idummy;
    if(verbose.debug)
      printf("DEBUG: Reading extended PuMa header version %d\n", version);
    if(version >= 1 && version <= 8) {
      idummy = (int *)dummy;
      dummy += sizeof(int);
      datafile->gentype = *idummy;
      if(verbose.debug) {
 printf("DEBUG: Extended header suggests that gentype=%d\n", datafile->gentype);
      }
      if(version >= 5) {
 bytevalue = *dummy;
 dummy += 1;
 datafile->xrangeset = bytevalue;
      }else {
 idummy = (int *)dummy;
 dummy += sizeof(int);
 datafile->xrangeset = *idummy;
      }
      fdummy = (float *)dummy;
      dummy += sizeof(float);
      datafile->xrange[0] = *fdummy;
      fdummy = (float *)dummy;
      dummy += sizeof(float);
      datafile->xrange[1] = *fdummy;
      if(version >= 5) {
 bytevalue = *dummy;
 dummy += 1;
 datafile->yrangeset = bytevalue;
      }else {
 idummy = (int *)dummy;
 dummy += sizeof(int);
 datafile->yrangeset = *idummy;
      }
      fdummy = (float *)dummy;
      dummy += sizeof(float);
      datafile->yrange[0] = *fdummy;
      fdummy = (float *)dummy;
      dummy += sizeof(float);
      datafile->yrange[1] = *fdummy;
      if(version >= 2) {
 if(version >= 5) {
   bytevalue = *dummy;
   dummy += 1;
   datafile->poltype = bytevalue;
 }else {
   idummy = (int *)dummy;
   dummy += sizeof(int);
   datafile->poltype = *idummy;
 }
 if(version >= 5) {
   bytevalue = *dummy;
   dummy += 1;
   datafile->feedtype = bytevalue;
 }else {
   idummy = (int *)dummy;
   dummy += sizeof(int);
   datafile->feedtype = *idummy;
 }
      }
      if(version >= 3) {
 if(version >= 5) {
   bytevalue = *dummy;
   dummy += 1;
   datafile->isDeFarad = bytevalue;
 }else {
   idummy = (int *)dummy;
   dummy += sizeof(int);
   datafile->isDeFarad = *idummy;
 }
      }
      if(version >= 4) {
 if(version >= 5) {
   bytevalue = *dummy;
   dummy += 1;
   datafile->isDePar = bytevalue;
 }else {
   idummy = (int *)dummy;
   dummy += sizeof(int);
   datafile->isDePar = *idummy;
 }
      }
      if(version >= 5) {
 bytevalue = *dummy;
 dummy += 1;
 datafile->cableSwap = bytevalue;
 bytevalue = *dummy;
 dummy += 1;
 datafile->cableSwapcor = bytevalue;
      }
      if(version >= 6) {
 fdummy = (float *)dummy;
 dummy += sizeof(float);
 datafile->telescope_X = *fdummy;
 fdummy = (float *)dummy;
 dummy += sizeof(float);
 datafile->telescope_Y = *fdummy;
 fdummy = (float *)dummy;
 dummy += sizeof(float);
 datafile->telescope_Z = *fdummy;
      }
      if(version >= 7) {
 fdummy = (float *)dummy;
 dummy += sizeof(float);
 datafile->freq_ref = *fdummy;
      }
      if(version >= 8) {
 bytevalue = *dummy;
 dummy += 1;
 datafile->isDebase = bytevalue;
      }
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readWSRTHeader: Undefined version of the extended header (%d).", version);
      return 0;
    }
  }
  return 1;
}
void FillPuMaHeader(Header_type *hdr, int obsID, int timefilenr, int freqband, int NFiles, int FileNum, int nrTimeSamples, int DataMJD, double DataTime, int StartMJD, char *ObsName, int StTime, double MaserOffset, double Dur, char *pulsar, double RA, double Dec, char *Epoch, Band_type *bands, int *BandsToClusMap, float *XPolScaleFac, int FIRFactor, int NSHARCsAdded, int NSampsAdded, int NFreqInFile, int Tsamp, int polmode, double AdjustInterval, int BitsPerSamp, double SigmaRange, int NrFrames, char *observatoryname, double longitude, double latitude, double height)
{
  int i;
  memset(hdr, 0, sizeof(Header_type));
  sprintf(hdr->gen.HdrVer, "DPC_1.3");
  sprintf(hdr->gen.Platform, "PuMaISim");
  sprintf(hdr->gen.ThisFileName, "%d.%05d.%d.puma", obsID, timefilenr, freqband);
  sprintf(hdr->gen.ScanNum, "%d", obsID);
  sprintf(hdr->gen.Comment, "PuMaISim");
  hdr->gen.NFiles=NFiles;
  hdr->gen.FileNum=FileNum;
  hdr->gen.TapeID[0] = 0;
  hdr->gen.NTapes = 0;
  hdr->gen.TapeNum = 0;
  hdr->gen.ParBlkSize = NrFrames*polmode*sizeof(Adjustments);
  hdr->gen.DataBlkSize = nrTimeSamples*BitsPerSamp*polmode*NFreqInFile/8;
  hdr->gen.Cluster[0] = FALSE;
  hdr->gen.Cluster[1] = FALSE;
  hdr->gen.Cluster[2] = FALSE;
  hdr->gen.Cluster[3] = FALSE;
  hdr->gen.Cluster[4] = FALSE;
  hdr->gen.Cluster[5] = FALSE;
  hdr->gen.Cluster[6] = FALSE;
  hdr->gen.Cluster[7] = FALSE;
  hdr->gen.Cluster[freqband] = TRUE;
  hdr->gen.DataMJD = DataMJD;
  hdr->gen.DataTime = DataTime;
  sprintf(hdr->obsy.Name, observatoryname);
  hdr->obsy.Long = longitude;
  hdr->obsy.Lat = latitude;
  hdr->obsy.Height = height;
  sprintf(hdr->obs.ObsName, ObsName);
  sprintf(hdr->obs.ObsType, "observation");
  sprintf(hdr->obs.LastCalID, "CalId");
  hdr->obs.StMJD = StartMJD;
  hdr->obs.StTime = StTime;
  hdr->obs.StLST = 0;
  hdr->obs.MaserOffset = MaserOffset;
  hdr->obs.Dur = Dur;
  hdr->obs.IonRM = 0;
  sprintf(hdr->src.Pulsar, "%s", pulsar);
  i = -1;
  do{
    i++;
    if(i == NAMELEN)
      break;
    if(hdr->src.Pulsar[i] == ' ')
      hdr->src.Pulsar[i] = '_';
  }while(hdr->src.Pulsar[i] != 0);
  hdr->src.RA = RA;
  hdr->src.Dec = Dec;
  sprintf(hdr->src.Epoch, "%s", Epoch);
  for(i = 0; i < MAXTELESCOPES; i++) {
    hdr->WSRT.Tel[i].Active = FALSE;
    hdr->WSRT.Tel[i].NoiseSrcOff = FALSE;
    hdr->WSRT.Tel[i].FrontendID[0] = 0;
    hdr->WSRT.Tel[i].FrontendStatus[0] = 0;
    hdr->WSRT.Tel[i].Tsys = 0;
  }
  for(i = 0; i < MAXFREQBANDS; i++) {
    memcpy(&(hdr->WSRT.Band[i].Operational), &(bands[i]), sizeof(Band_type));
  }
  hdr->WSRT.Backend[0] = 0;
  hdr->WSRT.ResFringeOff = TRUE;
  hdr->WSRT.AddBoxOn = TRUE;
  for(i = 0; i < MAXFREQBANDS; i++) {
    hdr->WSRT.BandsToClusMap[i] = BandsToClusMap[i];
  }
  hdr->mode.Nr = 1;
  for(i = 0; i < MAXFREQBANDS; i++) {
    hdr->mode.XPolScaleFac[i] = XPolScaleFac[i];
    hdr->mode.ActCluster[i] = bands[i].Operational;
  }
  hdr->mode.FIRFactor = FIRFactor;
  hdr->mode.NSHARCsAdded = NSHARCsAdded;
  hdr->mode.NSampsAdded = NSampsAdded;
  hdr->mode.NFreqInFile = NFreqInFile;
  hdr->mode.Tsamp = Tsamp;
  hdr->mode.Iout = FALSE;
  hdr->mode.Qout = FALSE;
  hdr->mode.Vout = FALSE;
  hdr->mode.Uout = FALSE;
  hdr->mode.Xout = FALSE;
  hdr->mode.Yout = FALSE;
  if(polmode == 1)
    hdr->mode.Iout = TRUE;
  else if(polmode == 2) {
    hdr->mode.Xout = TRUE;
    hdr->mode.Yout = TRUE;
  }else if(polmode == 4) {
    hdr->mode.Iout = TRUE;
    hdr->mode.Qout = TRUE;
    hdr->mode.Uout = TRUE;
    hdr->mode.Vout = TRUE;
  }
  hdr->mode.NDMs = 0;
  for(i = 0; i < MAXDMS; i++)
    hdr->mode.DM[i] = 0.000000;
  hdr->mode.RM=0;
  hdr->mode.DC_Dynamic = TRUE;
  hdr->mode.ScaleDynamic = TRUE;
  hdr->mode.AdjustInterval = AdjustInterval;
  hdr->mode.FloatsOut=FALSE;
  hdr->mode.BitsPerSamp = BitsPerSamp;
  hdr->mode.SigmaRange = SigmaRange;
  sprintf(hdr->software.Master, "Created by PuMaISim");
  for(i = 0; i < MAXDSPSPERBOARD; i++)
    hdr->software.DSP[i][0] = 0;
  sprintf(hdr->software.Filter, "Created by PuMaISim");
  hdr->check.ExitCodeObsy = 0;
  hdr->check.ExitCodePuMa = 0;
  for(i = 0; i < MAXFREQBANDS; i++)
    hdr->check.ExitCodeClust[i] = 0;
  hdr->check.ExitCodeDataConsistency = 0;
}
void beadj_convert_write(Adjustments outpadj, FILE *fout)
{
  int r;
  Adjustments beadj;
#ifdef __linux__
  PutBEint(&beadj.framenumber, outpadj.framenumber);
  PutBEfloat(&beadj.scale, outpadj.scale);
  PutBEfloat(&beadj.offset, outpadj.offset);
  r = fwrite(&beadj,sizeof(Adjustments),1,fout);
#endif
#ifdef __alpha
  PutBEint(&beadj.framenumber, outpadj.framenumber);
  PutBEfloat(&beadj.scale, outpadj.scale);
  PutBEfloat(&beadj.offset, outpadj.offset);
  r = fwrite(&beadj,sizeof(Adjustments),1,fout);
#endif
#ifdef __hpux
  r = fwrite(&outpadj,sizeof(Adjustments),1,fout);
#endif
}
void FillPuMaHeaderSimple(Header_type *hdr, int obsID, int timefilenr, int freqband, int flipped, float BW, float SkyMidFreq, long nrTimeSamples, double DataMJD, double StartMJD, char *ObsName, double Dur, char *pulsar, double RA, double Dec, int NSHARCsAdded, int NSampsAdded, int NFreqInFile, int Tsamp, int polmode, double AdjustInterval, int BitsPerSamp, long NrFrames, char *observatoryname, double longitude, double latitude, double height)
{
  int DataMJDint, StartMJDint, i, BandsToClusMap[8];
  float XPolScaleFac[8];
  char Epoch[NAMELEN];
  Band_type bands[8];
  Boolean NonFlip;
  DataMJDint = DataMJD;
  DataMJD -= DataMJDint;
  StartMJDint = StartMJD;
  StartMJD -= StartMJDint;
  StartMJD *= 24*3600;
  sprintf(Epoch, "J2006");
  if(flipped == 0)
    NonFlip = TRUE;
  else
    NonFlip = FALSE;
  for(i = 0; i < 8; i++) {
    bands[i].Operational = FALSE;
    bands[i].NonFlip = NonFlip;
    bands[i].SkyMidFreq = SkyMidFreq;
    bands[i].Width = BW;
    bands[i].MidFreq = 1.25;
  }
  bands[freqband].Operational = TRUE;
  BandsToClusMap[0] = 0;
  BandsToClusMap[1] = 1;
  BandsToClusMap[2] = 2;
  BandsToClusMap[3] = 3;
  BandsToClusMap[4] = 4;
  BandsToClusMap[5] = 5;
  BandsToClusMap[6] = 6;
  BandsToClusMap[7] = 7;
  XPolScaleFac[0] = 1;
  XPolScaleFac[1] = 1;
  XPolScaleFac[2] = 1;
  XPolScaleFac[3] = 1;
  XPolScaleFac[4] = 1;
  XPolScaleFac[5] = 1;
  XPolScaleFac[6] = 1;
  XPolScaleFac[7] = 1;
  FillPuMaHeader(hdr, obsID, timefilenr, freqband, 1, 1, nrTimeSamples, DataMJDint, DataMJD, StartMJDint, ObsName, StartMJD, 0, Dur, pulsar, RA, Dec, Epoch, bands, BandsToClusMap, XPolScaleFac, 0, NSHARCsAdded, NSampsAdded, NFreqInFile, Tsamp, polmode, AdjustInterval, BitsPerSamp, 3, NrFrames, observatoryname, longitude, latitude, height);
}
int writeWSRTHeader(datafile_definition datafile, verbose_definition verbose)
{
  Header_type puma_hdr;
  float *fdummy;
  int *idummy;
  char *dummy;
  if(datafile.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeWSRTHeader: Writing this data format is only implemented when the frequency channels are uniformely separated.");
    return 0;
  }
  FillPuMaHeaderSimple(&puma_hdr, 200700000, 1, 0, 0, get_bandwidth(datafile, verbose), get_centre_frequency(datafile, verbose), datafile.NrBins, datafile.mjd_start, datafile.mjd_start, "", get_tobs(datafile, verbose), datafile.psrname, datafile.ra, datafile.dec, 1, 1, datafile.NrFreqChan, get_tsamp(datafile, 0, verbose), datafile.NrPols, 0, 8, 1, datafile.observatory, observatory_long_geodetic(datafile), observatory_lat_geodetic(datafile), observatory_height_geodetic(datafile));
  puma_hdr.gen.ParBlkSize = 0;
  strncpy(puma_hdr.gen.ScanNum, datafile.scanID, NAMELEN);
  puma_hdr.redn.NBins = datafile.NrBins;
  double period;
  int ret;
  if(datafile.isFolded) {
    ret = get_period(datafile, 0, &period, verbose);
    if(ret == 2) {
      printerror(verbose.debug, "ERROR writeWSRTHeader (%s): Cannot obtain period", datafile.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = -1;
  }
  if(datafile.gentype == GENTYPE_SEARCHMODE || datafile.isFolded != 1 || ret == 1 || period < 0) {
    puma_hdr.redn.NTimeInts = datafile.NrBins*datafile.NrSubints;
    puma_hdr.redn.NBins = puma_hdr.redn.NTimeInts;
    puma_hdr.redn.GenType = 0;
    puma_hdr.redn.Folded = 0;
  }else {
    puma_hdr.redn.NTimeInts = datafile.NrSubints;
    puma_hdr.redn.GenType = 4;
    puma_hdr.redn.Folded = 1;
  }
  puma_hdr.redn.FoldPeriod = datafile.fixedPeriod;
  puma_hdr.redn.DeltaTime = get_tsamp(datafile, 0, verbose);
  puma_hdr.redn.MJDint = (int)datafile.mjd_start;
  puma_hdr.redn.MJDfrac = datafile.mjd_start - puma_hdr.redn.MJDint;
  puma_hdr.redn.FreqCent = get_centre_frequency(datafile, verbose);
  double chanbw;
  if(get_channelbandwidth(datafile, &chanbw, verbose) == 0) {
    printerror(verbose.debug, "ERROR writeWSRTHeader (%s): Cannot obtain channel bandwidth.", datafile.filename);
    return 0;
  }
  puma_hdr.redn.DeltaFreq = chanbw;
  puma_hdr.redn.DM = datafile.dm;
  puma_hdr.mode.RM = datafile.rm;
  puma_hdr.redn.NFreqs = datafile.NrFreqChan;
  puma_hdr.redn.Raw = FALSE;
  puma_hdr.redn.IsDedisp = datafile.isDeDisp;
  if(datafile.NrPols == 1) {
    puma_hdr.redn.OI = TRUE;
  }else if(datafile.NrPols == 2) {
    puma_hdr.redn.OX = TRUE;
    puma_hdr.redn.OY = TRUE;
  }else if(datafile.NrPols == 4) {
    puma_hdr.redn.OI = TRUE;
    puma_hdr.redn.OQ = TRUE;
    puma_hdr.redn.OU = TRUE;
    puma_hdr.redn.OV = TRUE;
  }else if(datafile.NrPols == 3) {
    puma_hdr.redn.OI = TRUE;
    puma_hdr.redn.OP = TRUE;
    puma_hdr.redn.OV = TRUE;
  }else if(datafile.NrPols == 5 && datafile.poltype == POLTYPE_ILVPAdPA) {
    puma_hdr.redn.OI = TRUE;
    puma_hdr.redn.OP = TRUE;
    puma_hdr.redn.OV = TRUE;
    puma_hdr.redn.OTheta = TRUE;
    puma_hdr.redn.Op = TRUE;
    puma_hdr.redn.OU = FALSE;
    puma_hdr.redn.OX = FALSE;
    puma_hdr.redn.OY = FALSE;
  }else if(datafile.NrPols == 8 && datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
    puma_hdr.redn.OI = TRUE;
    puma_hdr.redn.OP = TRUE;
    puma_hdr.redn.OV = TRUE;
    puma_hdr.redn.OTheta = TRUE;
    puma_hdr.redn.Op = TRUE;
    puma_hdr.redn.OU = TRUE;
    puma_hdr.redn.OX = TRUE;
    puma_hdr.redn.OY = TRUE;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeWSRTHeader: Cannot write out %ld polarization channels", datafile.NrPols);
    return 0;
  }
  dummy = puma_hdr.redn.ReductionDummy;
  strcpy(dummy, "psoft");
  dummy += 5;
  idummy = (int *)dummy;
  dummy += sizeof(int);
  *idummy = 8;
  idummy = (int *)dummy;
  dummy += sizeof(int);
  *idummy = datafile.gentype;
  *dummy = datafile.xrangeset;
  dummy += 1;
  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.xrange[0];
  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.xrange[1];
  *dummy = datafile.yrangeset;
  dummy += 1;
  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.yrange[0];
  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.yrange[1];
  *dummy = datafile.poltype;
  dummy += 1;
  *dummy = datafile.feedtype;
  dummy += 1;
  *dummy = datafile.isDeFarad;
  dummy += 1;
  *dummy = datafile.isDePar;
  dummy += 1;
  *dummy = datafile.cableSwap;
  dummy += 1;
  *dummy = datafile.cableSwapcor;
  dummy += 1;
  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.telescope_X;
  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.telescope_Y;
  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.telescope_Z;
  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.freq_ref;
  *dummy = datafile.isDebase;
  dummy += 1;
  rewind(datafile.fptr_hdr);
  pwheader(&puma_hdr, datafile.fptr_hdr);
  return 1;
}
long long PuMaFilepos(long long PulseNumber, long long FreqChan, long long polarization, long long binnr, long long NrSubints, long long NrBins, long long NrFreqChan, long long data_start)
{
  long long puma_filepos;
  puma_filepos = polarization*NrFreqChan*NrSubints*NrBins;
  puma_filepos += FreqChan*NrSubints*NrBins;
  puma_filepos += PulseNumber*NrBins;
  puma_filepos += binnr;
  puma_filepos *= sizeof(float);
  puma_filepos += data_start;
  return puma_filepos;
}
int readPulseWSRTData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse)
{
  long long filepos;
  filepos = polarization*datafile.NrFreqChan*datafile.NrSubints*datafile.NrBins;
  filepos += freq*datafile.NrSubints*datafile.NrBins;
  filepos += pulsenr*datafile.NrBins;
  filepos += binnr;
  filepos *= sizeof(float);
  filepos += datafile.datastart;
  fseeko(datafile.fptr, filepos, SEEK_SET);
  pumaread(pulse, sizeof(float), nrSamples, datafile.fptr);
  return 1;
}
int writePulseWSRTData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse)
{
  long long filepos;
  filepos = polarization*datafile.NrFreqChan*datafile.NrSubints*datafile.NrBins;
  filepos += freq*datafile.NrSubints*datafile.NrBins;
  filepos += pulsenr*datafile.NrBins;
  filepos += binnr;
  filepos *= sizeof(float);
  filepos += datafile.datastart;
  fseeko(datafile.fptr, filepos, SEEK_SET);
  pumawrite(pulse, sizeof(float), nrSamples, datafile.fptr);
  return 1;
}
int writePuMafile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, p;
  fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
  for(p = 0; p < datafile.NrPols; p++) {
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(n = 0; n < datafile.NrSubints; n++) {
 if(verbose.verbose && verbose.nocounters == 0) printf("writePuMafile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
 pumawrite(&data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], sizeof(float), datafile.NrBins, datafile.fptr);
      }
    }
  }
  if(verbose.verbose) printf("  Writing is done.              \n");
  return 1;
}
int readPuMafile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, p;
  if(verbose.verbose) {
    printf("Start reading PuMa file\n");
  }
  fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
  for(p = 0; p < datafile.NrPols; p++) {
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(n = 0; n < datafile.NrSubints; n++) {
 if(verbose.verbose && verbose.nocounters == 0)
   printf("  Progress reading PuMa file (%.1f%%)\r", 100.0*(n+(f+p*datafile.NrFreqChan)*datafile.NrSubints)/(float)(datafile.NrSubints*datafile.NrFreqChan*datafile.NrPols));
 pumaread(&data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], sizeof(float), datafile.NrBins, datafile.fptr);
      }
    }
  }
  if(verbose.verbose) printf("  Reading is done.                           \n");
  return 1;
}
int writeHistoryPuma(datafile_definition datafile, verbose_definition verbose)
{
  Header_type puma_hdr;
  char *txt, *txt_date, *user, *hostname, *txt_cmd;
  char questionmark[2];
  datafile_history_entry_definition *curHistoryEntry;
  questionmark[0] = '?';
  questionmark[1] = 0;
  txt = malloc(10000);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "Error writeHistoryPuma: Memory allocation error");
    return 0;
  }
  curHistoryEntry = &(datafile.history);
  while(curHistoryEntry->nextEntry != NULL) {
    curHistoryEntry = curHistoryEntry->nextEntry;
  }
  if(curHistoryEntry->timestamp != NULL)
    txt_date = curHistoryEntry->timestamp;
  else
    txt_date = questionmark;
  if(curHistoryEntry->user != NULL)
    user = curHistoryEntry->user;
  else
    user = questionmark;
  if(curHistoryEntry->hostname != NULL)
    hostname = curHistoryEntry->hostname;
  else
    hostname = questionmark;
  if(curHistoryEntry->cmd != NULL)
    txt_cmd = curHistoryEntry->cmd;
  else
    txt_cmd = questionmark;
  if(datafile.opened_flag) {
    sprintf(txt, "%s %s@%s: %s", txt_date, user, hostname, txt_cmd);
    memset(puma_hdr.redn.Command, 0, CMDLEN);
    strncpy(puma_hdr.redn.Command, txt, CMDLEN-1);
    fseeko(datafile.fptr_hdr, puma_hdr.redn.Command-(char *)(&puma_hdr), SEEK_SET);
    fwrite(puma_hdr.redn.Command, 1, CMDLEN, datafile.fptr_hdr);
    fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
    fseeko(datafile.fptr_hdr, datafile.datastart, SEEK_SET);
  }else {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING writeHistoryPuma: File not opened, ignoring command");
  }
  free(txt);
  return 1;
}
int readHistoryPuma(datafile_definition *datafile, verbose_definition verbose)
{
  Header_type puma_hdr;
  datafile_history_entry_definition *curHistoryEntry;
  if(datafile->opened_flag == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryPuma: File already closed?");
    return 0;
  }
  fseeko(datafile->fptr_hdr, 0, SEEK_SET);
  prheader(&puma_hdr, datafile->fptr_hdr);
  curHistoryEntry = &(datafile->history);
  curHistoryEntry->timestamp = NULL;
  curHistoryEntry->cmd = NULL;
  curHistoryEntry->user = NULL;
  curHistoryEntry->hostname = NULL;
  curHistoryEntry->nextEntry = NULL;
  curHistoryEntry->cmd = malloc(strlen(puma_hdr.redn.Command)+1);
  if(curHistoryEntry->cmd == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryPuma: Memory allocation error");
    return 0;
  }
  strcpy(curHistoryEntry->cmd, puma_hdr.redn.Command);
  return 1;
}
