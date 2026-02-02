C-------------------------------------------------------------------------------
C Spheral's initialization methods for the ANEOS EOS package.
C 
C Note I'm trying to reverse engineer ANEOS in all of this, so take it with a 
C caveat!
C
C Arguments:
C    filename : name of the ANEOS input file, normally something like ANEOS.INPUT
C    num      : number of materials 
C    izetl    : array of num EOS numbers

C-------------------------------------------------------------------------------
      subroutine ANEOS_initialize(in_filename, out_filename, num, izetl)

      implicit none

C ANEOS common blocks
      integer KLST, KINP
      COMMON /FILEOS/ KLST, KINP

C Declare the arguments.
      character*256 in_filename, out_filename, lineinp
      integer num, izetl(21)

C We'll use a common block to build a mapping from material numbers to offsets
C to help us in call_ANEOS.
      integer maxmat
      parameter (maxmat = 10)
      integer matnums(maxmat)
      common /Spheral_ANEOS_params/ matnums

C Local variables.
      integer i

C Build the offset material list, and check the input.
      do 10 i = 1, num
         if (izetl(i) .ge. 0) then
            print *, "ANEOS_initialize ERROR: make sure all EOS",
     &           " numbers are negative in initialize list."
            stop
         end if
         matnums(i) = -izetl(i)
 10   continue

C Open the files.
      open(10, file=in_filename, status='old')
      open(12, file=out_filename)
C      open(12, file=out_filename, status='new')
      kinp = 10
      klst = 12

C Call the main ANEOS initialization method.
      call ANEOS2(1, num, 0, izetl)

      end

C-------------------------------------------------------------------------------
C call_ANEOS
C
C Provides a wrapper around the ANEOS method used to compute the EOS response.
C-------------------------------------------------------------------------------
      subroutine call_ANEOS(matnum, T, rho, P, e, s, cv, dpdt, dpdr, cs)

      implicit none

      integer matnum
      double precision T, rho, P, e, s, cv, dpdt, dpdr, cs

c$$$      integer kpa
c$$$      double precision fkro
c$$$      call ANEOSD(T, rho, P, e, s, cv, dpdt, dpdr, fkro, cs, kpa, 
c$$$     &     matnum)
c$$$
c$$$      print *, " --> ", T, rho, P, e, s, cv, dpdt, dpdr, fkro, cs, kpa

C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      integer IPSQTS,MATBUF
      PARAMETER (MATBUF=64)
      real*8 sqts(MATBUF)
      COMMON /ANESQT/ SQTS,IPSQTS

C We'll use a common block to build a mapping from material numbers to offsets
C to help us in call_ANEOS.
      integer maxmat
      parameter (maxmat = 10)
      integer matnums(maxmat)
      common /Spheral_ANEOS_params/ matnums

C Local variables.
      integer i, matoffset

C Find the material number
      i = 0
 10   i = i + 1
      if (matnums(i) .eq. matnum) goto 20
      if (i .eq. maxmat) then
         print *, "call_ANEOS ERROR: unable to find material ", matnum
         stop
      end if
      goto 10

 20   ipsqts = 1
      sqts(ipsqts) = dsqrt(T)
      matoffset = 99*(i - 1)
      
      call ANEOS1(T, rho, P, e, s, cv, dpdt, dpdr, matoffset)

C ANEOS1 didn't compute the sound speed, so heres something.      
      cs = dsqrt(max(1e-10, dpdr))

      end

C-------------------------------------------------------------------------------
C get_ANEOS_atomicWeight
C
C Use the stored ANEOS data to compute the atomic weight of the given material.
C This has to be called after the initialization method above.
C-------------------------------------------------------------------------------
      function get_ANEOS_atomicWeight(matnum)
      
C Terrible idea!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision get_ANEOS_atomicWeight
      double precision minAtomicFrac
      integer matnum

C We'll use a common block to build a mapping from material numbers to offsets
C to help us in call_ANEOS.
      integer maxmat
      parameter (maxmat = 10)
      integer matnums(maxmat)
      common /Spheral_ANEOS_params/ matnums

C ANEOS stuff
      PARAMETER(NINPUT=48)      !NINPUT must be a multiple of 8!
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD

C Local variables.
      integer i, matoffset

C Find the material number
      i = 0
 110  i = i + 1
      if (matnums(i) .eq. matnum) goto 120
      if (i .eq. maxmat) then
         print *, "get_ANEOS_atomicWeight ERROR: unable to find ",
     &        " material ", matnum
         stop
      end if
      goto 110

 120  matoffset = 99*(i - 1)
      get_ANEOS_atomicWeight = ack(matoffset + 29)

C Find the minimum atomic fraction
      minAtomicFrac = 1.0
      matoffset = 30*(i - 1)
      do 150 i = 1, 30
         if (cot(i) .gt. 0.0) then
            minAtomicFrac = min(minAtomicFrac, cot(i))
         end if
 150  continue
      get_ANEOS_atomicWeight = get_ANEOS_atomicWeight/minAtomicFrac

      return
      end

C-------------------------------------------------------------------------------
C get_ANEOS_referenceDensity
C
C Use the stored ANEOS data to compute the reference mass density of the given
C  material.
C This has to be called after the initialization method above.
C-------------------------------------------------------------------------------
      function get_ANEOS_referenceDensity(matnum)
      
C Terrible idea!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision get_ANEOS_referenceDensity
      integer matnum

C We'll use a common block to build a mapping from material numbers to offsets
C to help us in call_ANEOS.
      integer maxmat
      parameter (maxmat = 10)
      integer matnums(maxmat)
      common /Spheral_ANEOS_params/ matnums

C ANEOS stuff
      PARAMETER(NINPUT=48)      !NINPUT must be a multiple of 8!
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD

C Local variables.
      integer i, matoffset

C Find the material number
      i = 0
 110  i = i + 1
      if (matnums(i) .eq. matnum) goto 120
      if (i .eq. maxmat) then
         print *, "get_ANEOS_atomicWeight ERROR: unable to find ",
     &        " material ", matnum
         stop
      end if
      goto 110

 120  matoffset = 99*(i - 1)
      get_ANEOS_referenceDensity = ack(matoffset + 11)
      return
      end

C-------------------------------------------------------------------------------
C get_ANEOS_atomicWeight
C
C Use the stored ANEOS data to compute the atomic weight of the given material.
C This has to be called after the initialization method above.
C
C I've cribbed and reimplemented the stuff in ANEOS2 for DIN(29) since it doesn't
C seem to be otherwise accessible.
C-------------------------------------------------------------------------------
c$$$      subroutine get_ANEOS_atomicWeight(matnum)
c$$$
c$$$      implicit none
c$$$
c$$$C Crap straight from ANEOS
c$$$      PARAMETER (MAXMAT=10)
c$$$      PARAMETER(NINPUT=48) !NINPUT must be a multiple of 8!
c$$$      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
c$$$     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
c$$$     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
c$$$     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
c$$$
c$$$      result = 0.0
c$$$      do 100 i = iz, izi
c$$$         IKK = ZZS(I)
c$$$         IKJ = IKK+(IKK*(IKK+1))/2
c$$$         TEIP = EIP(IKJ-IKK)           !Atomic weight of species I
c$$$         result = result + COT(I)*TEIP !mean atomic weight
c$$$ 100  continue
