!Colony batch run, allowing MPI parallel runs
! Compilation command: mpiifort ColonyBatchRun.f90 /O1 /Qopenmp
Module Variables
  IMPLICIT NONE
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER(I4B) :: ThreadID, NumThreads, NumFile 
  Character(Len = 99) :: ExecFile, ThreadIC
  Character(Len = 99), Allocatable :: InputFile(:)
End Module Variables
!
Program BatchRun
  USE Variables
  IMPLICIT NONE
!$  include 'mpif.h'
  INTEGER :: I, J, K
  INTEGER, Allocatable :: Idx(:)
  Character(Len = 99), Allocatable :: OutFileName(:)
  LOGICAL :: ExistQ
  EXTERNAL :: InitializeThreads, ReadFileName
  CALL InitializeThreads
  IF(ThreadID == 0) WRITE(*, '(A)')'Reading colony input file names...'
  CALL ReadFileName
  Allocate(Idx(NumFile), OutFileName(NumFile))
  K = 0
  DO I = 1, NumFile
    Inquire(FILE = Trim(InputFile(I)), Exist = ExistQ)
    IF(.NOT. ExistQ) CYCLE
    OPEN(10, File = Trim(InputFile(I)))
    READ(10, *) OutFileName(I)
    READ(10, *) OutFileName(I)
    OutFileName(I) = Trim(AdjustL(OutFileName(I)))
    CLOSE(10)
    Inquire(FILE = Trim(OutFileName(I)) // '.BestConfig', Exist = ExistQ)
    IF(ExistQ) CYCLE
    K = K + 1
    Idx(K) = I
  END DO
  IF(ThreadID == 0) WRITE(*, '(I4, A, I4, A, I4, A)') K, &
    ' files out of ', NumFile, ' to run by ', NumThreads, ' threads.'
!$ IF(NumThreads > 1) CALL MPI_BARRIER(MPI_COMM_WORLD, I)
  DO I = 1, K
    IF(MOD(I, NumThreads) /= ThreadID) CYCLE
    WRITE(*, '(A)')'Thread ' // Trim(ThreadIC) // &
      ' running data = ' // Trim(InputFile(Idx(I)))
    CALL SYSTEM(Trim(ExecFile) // ' SYS:1 IFN:' // Trim(InputFile(Idx(I))) &
      // ' TID:' // Trim(ThreadIC))
  END DO
!$ IF(NumThreads > 1) CALL MPI_BARRIER(MPI_COMM_WORLD, I)
!$ IF(NumThreads > 1) CALL MPI_FINALIZE(I)
End Program BatchRun
!
Subroutine InitializeThreads
  USE VARIABLES
  IMPLICIT NONE
!$  Include 'mpif.h'
!$  INTEGER :: I
  ThreadID = 0
  NumThreads = 1
!$  CALL MPI_INIT(I)
!$  IF(I /= MPI_SUCCESS)THEN
!$    WRITE(*,*)'Error in initialising MPI threads!'
!$    STOP
!$  ELSE
!$    CALL MPI_COMM_RANK(MPI_COMM_WORLD, ThreadID, I)
!$    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NumThreads, I)
!$  END IF
  WRITE(ThreadIC, '(I8)') ThreadID
  ThreadIC = Trim(AdjustL(ThreadIC))
END Subroutine InitializeThreads
!
SUBROUTINE ReadFileName
  USE Variables
  IMPLICIT NONE
  INTEGER(I4B) :: I
  CHARACTER(Len = 200) :: InFile, C200, Msg
  Logical :: ExistQ
  Msg = ''
  CALL Get_command_argument(1, ExecFile)
  IF(LEN_TRIM(ExecFile) <= 0) THEN
    Msg = 'No executable file name in command line!'
  ELSE
    CALL Get_command_argument(2, InFile)
    IF(LEN_TRIM(InFile) <= 0) Msg = &
      'No input file name in command line!'
  END IF
  IF(Len_Trim(Msg) == 0) THEN
    ExecFile = Trim(AdjustL(ExecFile))
    InFile = Trim(AdjustL(InFile))
    Inquire(FILE = Trim(ExecFile), Exist = ExistQ)
    IF(.NOT. ExistQ) THEN
      Msg = 'The executable, ' // Trim(ExecFile) // ', does not exist!'
    ELSE
      Inquire(FILE = Trim(InFile), Exist = ExistQ)
      IF(.NOT. ExistQ) Msg = 'The input file, ' // &
        Trim(InFile) // ', does not exist!'
    END IF
  END IF
  IF(Len_Trim(Msg) > 0) THEN
    IF(ThreadID == 0) WRITE(*, '(A)') Trim(Msg)
    IF(ThreadID == 0) WRITE(*, '(A)') 'Program stopping!'
    STOP
  END IF
  OPEN(10, File = TRIM(InFile))
  NumFile = 0
  DO
    READ(10, *, ERR=9, End=9) C200
    NumFile = NumFile + 1
  END DO
9 Close(10)
  IF(NumFile == 0) Msg = 'No colony input file names in file ' &
    // TRIM(InFile) // '!'
  Allocate(InputFile(NumFile))
  OPEN(10, File = TRIM(InFile))
  DO I = 1, NumFile
    READ(10, *) InputFile(I)
    InputFile(I) = Trim(AdjustL(InputFile(I)))
    Inquire(FILE = Trim(InputFile(I)), Exist = ExistQ)
    IF(ExistQ) CYCLE
    Msg = 'Colony input file, ' // Trim(InputFile(I)) // ', does not exist!'
    EXIT
  END DO
  CLOSE(10)
  IF(Len_Trim(Msg) > 0) THEN
    IF(ThreadID == 0) WRITE(*, '(A)') Trim(Msg)
    IF(ThreadID == 0) WRITE(*, '(A)') 'Program stopping!'
    STOP
  END IF
  CALL get_environment_variable("PATH", C200)
  IF(C200(1:1) == '/') ExecFile = './' // Trim(ExecFile)     !Linux/Apple
END SUBROUTINE ReadFileName
