SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4JSTEP(4)
C
C ----------------------------------------------------------------------C
C     Subroutine control
C ----------------------------------------------------------------------C
      INTEGER :: i, j, n_backstresses, it_num, converged, ind_alpha
C
C ----------------------------------------------------------------------C
C     Material properties
C ----------------------------------------------------------------------C
      REAL(8) :: elastic_modulus, Q_inf, b, D_inf, a,
     1shear_modulus, bulk_modulus, poission_ratio, mu2, lame_first
C
C ----------------------------------------------------------------------C
C     Used for intermediate calculations
C ----------------------------------------------------------------------C
      REAL(8) :: yield_stress, ep_eq, ep_eq_init, a_temp,
     1hard_iso_Q, hard_iso_D, hard_iso_total, a_dot_n,
     2plastic_mult, p_mult_numer, p_mult_denom, yield_function,
     3isotropic_modulus, kin_modulus,
     4stress_relative_norm, strain_trace, alpha_trace, e_k,
     5ID2_out_ID2, n_out_n, stress_hydro, sigma_vm,
     6Lam, n33_check, alpha_out_n, beta, theta_1, theta_2,
     7theta_3, srn2
C
C ----------------------------------------------------------------------C
C     Backstress arrays
C ----------------------------------------------------------------------C
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: alpha_k
      REAL(8), DIMENSION(:), ALLOCATABLE :: C_k, gamma_k
C
C ----------------------------------------------------------------------C
C     Tensors
C ----------------------------------------------------------------------C
      REAL(8), DIMENSION(6, 6) :: ID4, c_mat
      REAL(8), DIMENSION(6) :: strain_tens, strain_plastic,
     1yield_normal, alpha, strain_trial, stress_relative,
     2stress_dev, ID2, stress_tens, check, dstran_tens, alpha_diff,
     3alpha_upd
C
C ----------------------------------------------------------------------C
C     Parameters
C ----------------------------------------------------------------------C
      INTEGER :: N_BASIC_PROPS, TERM_PER_BACK, MAX_ITERATIONS,
     1I_ALPHA
      REAL(8) :: TOL, ONE, TWO, THREE, ZERO, SQRT23
C
C     N_BASIC_PROPS = 7
C     Material constants:
C       PROPS(1) = Elastic modulus
C       PROPS(2) = Poisson ratio
C       PROPS(3) = Initial yield stress
C       PROPS(4) = Q_inf
C       PROPS(5) = b
C       PROPS(6) = D_inf
C       PROPS(7) = a
C
C     Backstress constants:
C       PROPS(8)  = C_1
C       PROPS(9)  = gamma_1
C       PROPS(10) = C_2
C       PROPS(11) = gamma_2
C       ...
C
C     STATEV layout after modification:
C       STATEV(1)    = Equivalent plastic strain
C       STATEV(2:7)  = Internal plastic strain vector
C                      11, 22, 33, 12, 13, 23
C                      Shear components are engineering shear strains.
C       STATEV(8)    = True tensor plastic strain EP11
C       STATEV(9)    = True tensor plastic strain EP12
C       STATEV(10)   = True tensor plastic strain EP13
C       STATEV(11)   = True tensor plastic strain EP21
C       STATEV(12)   = True tensor plastic strain EP22
C       STATEV(13)   = True tensor plastic strain EP23
C       STATEV(14)   = True tensor plastic strain EP31
C       STATEV(15)   = True tensor plastic strain EP32
C       STATEV(16)   = True tensor plastic strain EP33
C       STATEV(17:)  = Backstress components
C
C     Therefore:
C       NSTATV = 16 + 6 * n_backstresses
C
      PARAMETER(TOL=1.0D-10,
     1N_BASIC_PROPS=7, TERM_PER_BACK=2, MAX_ITERATIONS=1000,
     2ONE=1.0D0, TWO=2.0D0, THREE=3.0D0, ZERO=0.D0,
     3SQRT23=SQRT(2.0D0/3.0D0), I_ALPHA=16)
C
C ----------------------------------------------------------------------C
C     Subroutine start
C ----------------------------------------------------------------------C
C
C     This UMAT assumes a full 3D stress-strain state with NTENS = 6.
C     The output true plastic strain tensor STATEV(8:16) is constructed
C     from the 6-component internal plastic strain vector.
C
      IF (NTENS .NE. 6) THEN
        PRINT *, "ERROR: This modified UMAT requires NTENS = 6."
        PRINT *, "Current NTENS = ", NTENS
        PRINT *, "Element = ", NOEL, ", Integration point = ", NPT
        CALL XIT
      END IF
C
C ----------------------------------------------------------------------C
C     Get the number of backstresses
C ----------------------------------------------------------------------C
      n_backstresses = (nprops - N_BASIC_PROPS) / TERM_PER_BACK
      IF (n_backstresses .EQ. 0) THEN
        PRINT *, "No backstresses defined, exiting!"
        CALL XIT
      END IF
C
C ----------------------------------------------------------------------C
C     Check number of state variables
C ----------------------------------------------------------------------C
      IF (NSTATV .LT. I_ALPHA + n_backstresses * NTENS) THEN
        PRINT *, "ERROR: Not enough state variables are defined."
        PRINT *, "Required NSTATV = ", I_ALPHA + n_backstresses * NTENS
        PRINT *, "Current  NSTATV = ", NSTATV
        PRINT *, "Please update *DEPVAR in the input file."
        CALL XIT
      END IF
C
C ----------------------------------------------------------------------C
C     Allocate the backstress related arrays
C ----------------------------------------------------------------------C
      ALLOCATE(C_k(n_backstresses))
      ALLOCATE(gamma_k(n_backstresses))
      ALLOCATE(alpha_k(n_backstresses, ntens))
C
C ----------------------------------------------------------------------C
C     Initialize
C ----------------------------------------------------------------------C
      ddsdde(:, :) = ZERO
      ID4(:, :) = ZERO
C
C     4th order symmetric identity tensor in Voigt form
      DO i = 1, ndi
        ID4(i, i) = ONE
      END DO
      DO i = ndi+1, ntens
        ID4(i, i) = ONE / TWO
      END DO
C
C     2nd order symmetric identity tensor
      ID2(:) = (/ ONE, ONE, ONE, ZERO, ZERO, ZERO /)
C
C ----------------------------------------------------------------------C
C     Read in state variables
C ----------------------------------------------------------------------C
C     STATEV(1)   = Equivalent plastic strain
C     STATEV(2:7) = Internal plastic strain vector
C                  11, 22, 33, 12, 13, 23
C                  Shear components are engineering shear strains.
C
C     STATEV(8:16) are only output variables for the full true plastic
C     strain tensor. They are reconstructed later from STATEV(2:7).
C
C     STATEV(17:) = Backstress components
C
      ep_eq = statev(1)
      ep_eq_init = statev(1)
C
      CALL ROTSIG(statev(2), drot, strain_plastic, 2, ndi, nshr)
C
      alpha(:) = ZERO
      DO i = 1, n_backstresses
        ind_alpha = I_ALPHA + 1 + (i - 1) * ntens
        CALL ROTSIG(statev(ind_alpha), drot, alpha_k(i, :),
     1  1, ndi, nshr)
        alpha = alpha + alpha_k(i, :)
      END DO
C
C ----------------------------------------------------------------------C
C     Read in the material properties
C ----------------------------------------------------------------------C
      elastic_modulus = props(1)
      poission_ratio = props(2)
      yield_stress = props(3)
      q_inf = props(4)
      b = props(5)
      d_inf = props(6)
      a = props(7)
C
      DO i = 1, n_backstresses
        c_k(i) = props((N_BASIC_PROPS - 1) + 2 * i)
        gamma_k(i) = props(N_BASIC_PROPS + 2 * i)
      END DO
C
C ----------------------------------------------------------------------C
C     Calculate elastic parameters
C ----------------------------------------------------------------------C
      shear_modulus = elastic_modulus / (TWO * (ONE + poission_ratio))
      bulk_modulus = elastic_modulus /
     1(THREE * (ONE - TWO * poission_ratio))
      mu2 = TWO * shear_modulus
C
C ----------------------------------------------------------------------C
C     Set-up strain tensor
C ----------------------------------------------------------------------C
C     Tensors are stored as:
C       1 = 11
C       2 = 22
C       3 = 33
C       4 = 12
C       5 = 13
C       6 = 23
C
      strain_tens = stran + dstran
C
C ----------------------------------------------------------------------C
C     Elastic trial step
C ----------------------------------------------------------------------C
C
C     Tensor of elastic moduli
      DO j = 1, ntens
        DO i = 1, ntens
          ID2_out_ID2 = ID2(i) * ID2(j)
          c_mat(i, j) = ID2_out_ID2 * bulk_modulus +
     1    mu2 * (ID4(i, j) - ONE / THREE * ID2_out_ID2)
        END DO
      END DO
C
C     Trial stress tensor
      stress_tens = MATMUL(c_mat, (strain_tens - strain_plastic))
C
      stress_hydro = SUM(stress_tens(1:3)) / THREE
      strain_trace = SUM(strain_tens(1:3))
C
      DO i = 1, ndi
        stress_dev(i) = stress_tens(i) - stress_hydro
        stress_relative(i) = stress_dev(i) - alpha(i)
      END DO
C
      DO i = ndi+1, ntens
        stress_dev(i) = stress_tens(i)
        stress_relative(i) = stress_dev(i) - alpha(i)
      END DO
C
      stress_relative_norm =
     1SQRT(dotprod6(stress_relative, stress_relative))
C
C ----------------------------------------------------------------------C
C     Yield condition
C ----------------------------------------------------------------------C
      hard_iso_Q = q_inf * (ONE - EXP(-b * ep_eq))
      hard_iso_D = d_inf * (ONE - EXP(-a * ep_eq))
      hard_iso_total = yield_stress + hard_iso_Q - hard_iso_D
C
      yield_function = stress_relative_norm - SQRT23 * hard_iso_total
C
      IF (yield_function .GT. TOL) THEN
        converged = 0
      ELSE
        converged = 1
      END IF
C
C     Calculate the normal to the yield surface
      yield_normal = stress_relative / (TOL + stress_relative_norm)
C
C ----------------------------------------------------------------------C
C     Radial return mapping if plastic loading
C ----------------------------------------------------------------------C
C
C     Calculate the consistency parameter, namely plastic multiplier
      plastic_mult = ZERO
      it_num = 0
C
      DO WHILE ((converged .EQ. 0) .AND.
     1(it_num .LT. MAX_ITERATIONS))
C
        it_num = it_num + 1
C
C       Calculate the isotropic hardening parameters
        hard_iso_Q = q_inf * (ONE - EXP(-b * ep_eq))
        hard_iso_D = d_inf * (ONE - EXP(-a * ep_eq))
        hard_iso_total = yield_stress + hard_iso_Q - hard_iso_D
C
        isotropic_modulus = b * (q_inf - hard_iso_Q) -
     1  a * (d_inf - hard_iso_D)
C
C       Calculate the kinematic hardening parameters
        kin_modulus = ZERO
        DO i = 1, n_backstresses
          e_k = EXP(-gamma_k(i) * (ep_eq - ep_eq_init))
          kin_modulus = kin_modulus + C_k(i) * e_k
     1    - SQRT(THREE / TWO) * gamma_k(i) * e_k
     2    * dotprod6(yield_normal, alpha_k(i, :))
        END DO
C
        a_dot_n = ZERO
        alpha_upd(:) = ZERO
C
        DO i = 1, n_backstresses
          e_k = EXP(-gamma_k(i) * (ep_eq - ep_eq_init))
          alpha_upd = alpha_upd + e_k * alpha_k(i, :)
     1    + SQRT23 * C_k(i) / gamma_k(i) * (ONE - e_k)
     2    * yield_normal
        END DO
C
C       n : delta_alpha
        a_dot_n = dotprod6(alpha_upd - alpha, yield_normal)
C
        p_mult_numer = stress_relative_norm -
     1  (a_dot_n + SQRT23 * hard_iso_total + mu2 * plastic_mult)
C
        p_mult_denom = -mu2 *
     1  (ONE + (kin_modulus + isotropic_modulus) /
     2  (THREE * shear_modulus))
C
C       Update plastic multiplier
        plastic_mult = plastic_mult - p_mult_numer / p_mult_denom
C
C       Update equivalent plastic strain
        ep_eq = ep_eq_init + SQRT23 * plastic_mult
C
        IF (ABS(p_mult_numer) .LT. TOL) THEN
          converged = 1
        END IF
C
      END DO
C
C ----------------------------------------------------------------------C
C     Update variables
C ----------------------------------------------------------------------C
      IF (it_num .EQ. 0) THEN
C
C       Elastic loading
        stress = stress_tens
C
      ELSE
C
C       Plastic loading
C
C       The normal components are updated directly.
C       The shear components in the internal strain vector are
C       engineering shear components, so they need the second update.
C
        strain_plastic = strain_plastic + plastic_mult * yield_normal
C
        strain_plastic(4:6) = strain_plastic(4:6)
     1  + plastic_mult * yield_normal(4:6)
C
        stress = MATMUL(c_mat, (strain_tens - strain_plastic))
C
        alpha_diff = alpha
        alpha(:) = ZERO
C
C       Update backstress components
        DO i = 1, n_backstresses
          e_k = EXP(-gamma_k(i) * (ep_eq - ep_eq_init))
          alpha_k(i, :) = e_k * alpha_k(i, :) +
     1    SQRT23 * yield_normal * C_k(i) / gamma_k(i)
     2    * (ONE - e_k)
          alpha = alpha + alpha_k(i, :)
        END DO
C
        alpha_diff = alpha - alpha_diff
C
      END IF
C
C ----------------------------------------------------------------------C
C     Tangent modulus
C ----------------------------------------------------------------------C
      IF (it_num .EQ. 0) THEN
C
C       Elastic loading
        DO j = 1, ntens
          DO i = 1, ntens
            ddsdde(i, j) = c_mat(i, j)
          END DO
        END DO
C
        DO j = ndi+1, ntens
          ddsdde(j, j) = shear_modulus
        END DO
C
      ELSE
C
C       Plastic loading
        beta = ONE +
     1  (kin_modulus + isotropic_modulus) /
     2  (THREE * shear_modulus)
C
        theta_1 = ONE - mu2 * plastic_mult / stress_relative_norm
        theta_3 = ONE / (beta * stress_relative_norm)
        theta_2 = ONE / beta
     1  + dotprod6(yield_normal, alpha_diff) * theta_3
     2  - (ONE - theta_1)
C
        DO j = 1, ntens
          DO i = 1, ntens
            ID2_out_ID2 = ID2(i) * ID2(j)
            n_out_n = yield_normal(i) * yield_normal(j)
            alpha_out_n = alpha_diff(i) * yield_normal(j)
C
            ddsdde(i, j) = bulk_modulus * ID2_out_ID2
     1      + mu2 * theta_1 *
     2      (ID4(i, j) - ONE / THREE * ID2_out_ID2)
     3      - mu2 * theta_2 * n_out_n
     4      + mu2 * theta_3 * alpha_out_n
C
          END DO
        END DO
C
        ddsdde = ONE / TWO * (TRANSPOSE(ddsdde) + ddsdde)
C
      END IF
C
C ----------------------------------------------------------------------C
C     Update the state variables
C ----------------------------------------------------------------------C
C
C     Equivalent plastic strain
      statev(1) = ep_eq
C
C     Internal plastic strain vector.
C     Components:
C       STATEV(2) = EP11
C       STATEV(3) = EP22
C       STATEV(4) = EP33
C       STATEV(5) = GAMMA_P12 = 2 * EP12
C       STATEV(6) = GAMMA_P13 = 2 * EP13
C       STATEV(7) = GAMMA_P23 = 2 * EP23
C
      DO i = 1, ntens
        statev(i + 1) = strain_plastic(i)
      END DO
C
C ----------------------------------------------------------------------C
C     New output: true full plastic strain tensor
C ----------------------------------------------------------------------C
C
C     The following STATEV variables are for output and post-processing.
C     They form the true tensor:
C
C       [ EP11  EP12  EP13 ]
C       [ EP21  EP22  EP23 ]
C       [ EP31  EP32  EP33 ]
C
C     Since the internal shear components are engineering shear strains:
C
C       EP12 = 0.5 * GAMMA_P12
C       EP13 = 0.5 * GAMMA_P13
C       EP23 = 0.5 * GAMMA_P23
C
      statev(8)  = strain_plastic(1)
      statev(9)  = HALF(strain_plastic(4))
      statev(10) = HALF(strain_plastic(5))
C
      statev(11) = HALF(strain_plastic(4))
      statev(12) = strain_plastic(2)
      statev(13) = HALF(strain_plastic(6))
C
      statev(14) = HALF(strain_plastic(5))
      statev(15) = HALF(strain_plastic(6))
      statev(16) = strain_plastic(3)
C
C ----------------------------------------------------------------------C
C     Update backstress state variables
C ----------------------------------------------------------------------C
      DO i = 1, n_backstresses
        DO j = 1, ntens
          statev(I_ALPHA + j + (i - 1) * ntens) = alpha_k(i, j)
        END DO
      END DO
C
C ----------------------------------------------------------------------C
C     Reduce time increment if did not converge
C ----------------------------------------------------------------------C
      IF (it_num .EQ. MAX_ITERATIONS) THEN
        PRINT *, "WARNING: Return mapping in integration point ", npt,
     1  " of element ", noel, " did not converge."
        PRINT *, "Reducing time increment to 1/10 of current value."
        PNEWDT = 0.10
      END IF
C
C ----------------------------------------------------------------------C
C     Deallocate arrays
C ----------------------------------------------------------------------C
      IF (ALLOCATED(C_k)) DEALLOCATE(C_k)
      IF (ALLOCATED(gamma_k)) DEALLOCATE(gamma_k)
      IF (ALLOCATED(alpha_k)) DEALLOCATE(alpha_k)
C
      RETURN
C
      CONTAINS
C
C ----------------------------------------------------------------------C
C     Define dot product for vectors
C ----------------------------------------------------------------------C
      PURE FUNCTION dotprod6(A, B) RESULT(C)
C
C     Returns the dot product of two symmetric length-6 vectors.
C     For tensor inner product:
C
C       A : B = A11 B11 + A22 B22 + A33 B33
C             + 2 A12 B12 + 2 A13 B13 + 2 A23 B23
C
      REAL(8), INTENT(IN) :: A(6), B(6)
      REAL(8)             :: C
      INTEGER             :: i
C
      C = 0.0D0
C
      DO i = 1, 3
        C = C + A(i) * B(i)
      END DO
C
      DO i = 4, 6
        C = C + TWO * A(i) * B(i)
      END DO
C
      END FUNCTION dotprod6
C
C ----------------------------------------------------------------------C
C     Half function
C ----------------------------------------------------------------------C
      PURE FUNCTION HALF(A) RESULT(B)
C
      REAL(8), INTENT(IN) :: A
      REAL(8)             :: B
C
      B = 0.5D0 * A
C
      END FUNCTION HALF
C
      END
