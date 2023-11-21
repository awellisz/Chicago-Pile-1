program fission
    implicit none 
    ! Number of neutrons to simulate
    integer :: num_neutrons
    ! Thickness of graphite slab, in cm
    integer :: slab_thickness 
    ! Counters
    integer :: num_absorbed, num_backscattered, num_thermalized, num_traversed
    ! File I/O
    integer :: io
    
    num_neutrons = 5000

    open(unit=io, file="test1.csv", status='replace', action='write')
    write(io,*)"slab_thickness,n_absorbed,n_backscattered,n_traversed,n_thermalized"

    do slab_thickness = 1, 51, 3
        print *,"Simulating D = ",slab_thickness,"..."
        call run_simulation(num_neutrons, slab_thickness, num_absorbed, num_backscattered, num_thermalized, num_traversed)
        ! Write to file 
        write(io, *)slab_thickness,",",num_absorbed,",",num_backscattered,",",num_traversed,",",num_thermalized
    end do 

end program fission

! energy_loss: Calculate the energy loss in collision, use isotropic 
!              (s-wave) approximation
!   energy_in:  input energy of particle
!   angle:      scattering angle
!   energy_out: final energy
subroutine energy_loss(energy_in, angle)
    implicit none 
    real, intent(inout) :: energy_in
    real, intent(out) :: angle 
    real :: energy_out
    real :: rn1
    real :: max_recoil_energy 
    real :: Z ! Atomic number
    real :: eta

    Z = 12 ! Graphite (C12)

    ! Approximation: round any energy below 0.1 eV to 0.02 eV
    if (energy_in < 1.e-7) then 
        energy_out = 2.e-8
    else
        ! Maximum recoil energy that can be imparted
        max_recoil_energy = ((4 * Z) / ((1 + Z)**2)) * energy_in
        call random_number(rn1)
        ! Probability distribution for recoil energies is flat from 0 to max
        max_recoil_energy = max_recoil_energy * rn1 
        energy_out = energy_in - max_recoil_energy
    end if

    if (energy_in < 1.e-7) then 
        call random_number(rn1)
        angle = acos((2. * rn1) - 1.)
    else 
        eta = acos(sqrt((max_recoil_energy/energy_in)*(((1.+Z)**2)/(4.*Z))))
        angle = atan((sin(2.*eta))/((1./Z) - cos(2.*eta)))
    end if

    ! Update energy
    energy_in = energy_out 
end subroutine energy_loss

! init_energy: Uses von Neumann's method to calculate the energy of a neutron 
!              emitted from a fission source
!   energy: the energy of the neutron (in MeV)
subroutine init_energy(energy)
    implicit none
    real, intent(out) :: energy
    real :: rn1, rn2
    real :: prob 

    prob = 0.0
    rn2 = 1.0

    do while (rn2 > prob)

        call random_number(rn1)
        ! P = 0 above ~10 MeV
        rn1 = rn1 * 9.999
        ! Apply Maxwell spectrum
        prob = sqrt(rn1) * exp(-rn1 / 1.4)
        
        call random_number(rn2)
        ! Maximum of the distribution as written is 0.5
        rn2 = rn2 * 0.5
    enddo

    energy = rn1 
end subroutine init_energy

! sample_watt: Samples the Watt distribution for Uranium-235 to choose a random
!              initial energy of a simulated neutron.
!   energy: the energy of the neutron (result in MeV)
! subroutine sample_watt(energy)
!     implicit none 
!     real, intent(out) :: energy 
!     real :: K, L, M, x, y, rn1, rn2 
!     real :: a, b 

!     ! Values for Uranium-235
!     a = 0.75
!     b = 1.0

!     K = 1 + (b / 8 * a)
!     L = (K + sqrt(K**2 - 1)) / a 
!     M = a*L - 1
!     do while(.true.)
!         call random_number(rn1)
!         call random_number(rn2)
!         x = -log(rn1)
!         y = -log(rn2)
!         if ((y - M*(x + 1))**2 < (b * L * x)) then 
!         E = L * x
!         exit
!         endif
!     enddo
! end subroutine sample_watt

! cross_section: Calculates the approximate cross section (in barns) for 
!                absorption and elastic scattering in carbon
! energy: energy in MeV
! xelastic: cross section for elastic scattering
! xabsorbed: cross section for absorption
subroutine cross_section(energy, xelastic, xabsorbed)
    real, intent(in) :: energy
    real, intent(out) :: xelastic
    real, intent(out) :: xabsorbed
    real :: energy_eV 

    ! We need energy in eV, not MeV
    eev = energy * 1.e6
    
    ! Approximations...
    if (eev < 1.e4) then 
        xelastic = 5.
    else 
        xelastic = 10.5 - (1.346 * log10(eev))
    end if 

    ! Convert to barns (1 barn is 10^{-24} cm^2)
    xelastic = xelastic * 1.e-24

    if (eev < 1.e3) then 
        xabsorbed = (6.442e-4) * (eev**(-0.49421))
    else if (eev < 1.e5) then 
        xabsorbed = 1.5e-5
    else if (eev < 5.e6) then 
        xabsorbed = 2.e-5 
    else 
        xabsorbed = (4.e-6) * exp(eev * 3.2189e-7)
    endif 

    xabsorbed = xabsorbed * 1.e-24

end subroutine cross_section

! euler: Takes the original linear trajectory, rotates it to lie along the z-axis,
!        generates a vector at zenith angle theta (= scattering angle) and 
!        azimuthal angle (= random * 2pi). The original axis is now rotated back
!        taking the scattering vector with it. Now we have the scattered direction
!        vector (Sx, Sy, Sz). We use Euler angles to perform the transformation.
!subroutine euler(Ex, Ey, Ez, angle, Sx, Sy, Sz)
subroutine euler(in_dir, angle)
    implicit none
    real, intent(inout) :: in_dir(3)
    !real, intent(inout) :: Ex, Ey, Ez 
    real, intent(in) :: angle 
    !real, intent(out) :: Sx, Sy, Sz
    real :: scatter_dir(3)

    ! real :: scatter_dir

    real :: Ex, Ey, Ez

    real :: s, arg, rn1
    real :: alpha, beta, gamma, theta, phi
    real :: sco1, sco2
    real :: r11, r12, r13, r21, r22, r23, r31, r32, r33 ! rotation matrix elements
    real :: S0x, S0y, S0z ! Pre-transformation propagation vector?

    ! Normalize direction to a unit vector (in case it wasn't)
    s = norm2(in_dir)
    in_dir = in_dir / s
    Ex = in_dir(1)
    Ey = in_dir(2)
    Ez = in_dir(3)

    beta = acos(Ez)

    ! Watch out for divide by zero errors here
    arg = Ey / sin(beta)
    if (abs(arg) > 1.0) then  ! huh??
        arg = arg / (1.0001 * abs(arg))
    end if 

    alpha = asin(arg)

    ! Approximation?
    if (abs(beta) < 0.027) then 
        alpha = 0.0
    end if 

    sco1 = abs(cos(alpha) * sin(beta) + Ex)
    sco2 = abs(Ex) 

    if (sco1 < sco2) then 
        beta = -beta 
        alpha = -alpha
    end if 

    gamma = 0.0 

    ! We now have the Euler angles of rotation between the z-axis to the
    !  direction of the initial particle
    theta = angle 
    call random_number(rn1)
    phi = 6.2831853 * rn1 

    ! We have now scattered the particle from the z-axis and must rotate it
    !  to the original unscattered particle direction. 
    ! Calculates rotation matrix
    r11 = cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma)
    r12 = cos(beta)*sin(alpha)*cos(gamma)+cos(alpha)*sin(gamma)
    r13 =-sin(beta)*cos(gamma)
    r21 =-sin(gamma)*cos(beta)*cos(alpha)-sin(alpha)*cos(gamma)
    r22 =-sin(gamma)*cos(beta)*sin(alpha)+cos(alpha)*cos(gamma)
    r23 = sin(beta)*sin(gamma)
    r31 = sin(beta)*cos(alpha)
    r32 = sin(alpha)*sin(beta)
    r33 = cos(beta)
    S0x = sin(theta)*cos(phi)
    S0y = sin(theta)*sin(phi)
    S0z = cos(theta)
    ! Unit propagation vector of the scattered particle in the original frame
    scatter_dir = [r11*S0x+r21*S0y+r31*S0z, r12*S0x+r22*S0y+r32*S0z, r13*S0x+r23*S0y+r33*S0z]

    in_dir = scatter_dir

end subroutine euler


subroutine run_simulation(num_neutrons, slab_thickness,  num_absorbed, &
                            num_backscattered, num_thermalized, num_traversed)
    implicit none
    integer, intent(in) :: num_neutrons, slab_thickness
    integer, intent(out) :: num_absorbed, num_backscattered, num_thermalized, num_traversed

    ! CONSTANTS 
    real :: NA   ! Avogadro's number
    real :: rho  ! Density of medium (graphite) 
    integer :: Z ! atomic number (carbon-12)
    real :: pi 

    ! PROPERTIES OF NEUTRONS
    real :: position(3)          ! [x, y, z] position of neutron 
    real :: theta, phi           ! Angles for generating velocity vector
    real :: velocity(3)          ! Unit velocity vector
    real :: scatter_angle        ! Scattering angle
    real :: scatter_dir          ! Scattering direction unit vector
    real :: energy               ! Energy
    real :: energy_out           ! Energy after a scattering event
    real :: xelastic, xabsorbed  ! Cross sections
    real :: dnext                ! Distance to the next interaction point

    ! MISC
    integer :: i, j  ! For iterating
    real :: rn       ! For random numbers
    
    pi = 4 * atan(1.0)
    NA = 6.023e23 
    rho = 2.16 ! g/cm^3
    Z = 12     ! g/mol

    !print *, "Hello world!"

    num_backscattered = 0
    num_traversed = 0
    num_thermalized = 0
    num_absorbed = 0

    ! Main simulation loop to do for each neutron
    simulation_loop: do i = 1, num_neutrons
        !print *,"Simulating neutron ",i
        call random_seed() 

        ! Initialize position to origin 
        position(:) = 0

        ! Initialize energy 
        call init_energy(energy)

        ! Build unit velocity vector
        call random_number(phi)
        phi = 2*pi * phi
        theta = -1.
        ! Make sure theta is positive
        do while (theta < 0)
            call random_number(theta)
            theta = acos(2*theta - 1) ! Proper isotropic direction
        ! Convert angles to cartesian vector
        velocity = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]

        ! Individual neutron loop
        ! (Exits loop if neutron backscatters/traverses/is absorbed)
        neutron_loop: do while (.true.)

            ! Find cross sections for elastic scattering / absorption
            call cross_section(energy, xelastic, xabsorbed)
            ! Determine distance to next interaction
            call random_number(rn)
            ! Depends on total cross section
            dnext = -log(rn) * Z / ((xelastic + xabsorbed) * rho * NA)

            ! Calculate new position (at interaction)
            position = position + dnext * velocity

            ! Check for backscattering (Z < 0)
            if (position(3) < 0.0) then 
                num_backscattered = num_backscattered + 1
                !print *,"Backscattered"
                cycle simulation_loop
            end if

            ! Check for traversal (Z > D)
            if (position(3) > slab_thickness) then 
                num_traversed = num_traversed + 1
                !print *,"Traversed"
                ! Was it thermalized? (energy below 500 eV)
                if ((energy * 1.e6) < 500) then
                    !print *,"Thermalized"
                    num_thermalized = num_thermalized + 1
                end if 
                cycle simulation_loop
            end if

            ! Check for absorption 
            call random_number(rn)
            ! Calculate robability of absorption as ratio of cross sections
            if (rn < (xabsorbed / (xabsorbed + xelastic))) then 
                !print *,"Absorbed"
                num_absorbed = num_absorbed + 1
                cycle simulation_loop
            end if

            ! If not backscattered/traversed/absorbed, calculate new energy
            !  after energy loss and scattering angle.
            call energy_loss(energy, scatter_angle)

            ! Find the new trajectory
            call euler(velocity, scatter_angle)

        end do neutron_loop
    end do simulation_loop

end subroutine run_simulation