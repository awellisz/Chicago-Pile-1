program pile
    implicit none 
    ! Number of neutrons to simulate
    integer :: num_neutrons
    ! Characteristic thickness
    real :: D
    ! Neutron energy (for simulate_energy)
    real :: energy
    ! Counters
    real :: total_fission_xsects
    integer :: num_absorbed
    ! misc
    integer :: i, j

    ! u235fission.in cross section vs. neutron energy data
    real :: xsections(31882)
    real :: energies(31882)


    ! Read data from u235fission.in
    open(unit=1, file="u235fission.in", status="old", action="read")
    do i = 1, 31882
      read(1, *) energies(i), xsections(i)
    end do
    close(1)

    open(unit=2, file="pile_xsections.csv", status="replace", action="write")
    write(2, *)"D,fission_xsection,num_absorbed"

    num_neutrons = 100000
    total_fission_xsects = 0.0
    num_absorbed = 0
    ! D = 10. ! cm

    do i = 1, 70
        D = 1. * i
        print *, "Simulating D = ",D
        call run_simulation(num_neutrons, D, num_absorbed, total_fission_xsects, &
                          energies, xsections)
        print *,"Fission cross sections:",total_fission_xsects
        write(2,*)D,",",total_fission_xsects,",",num_absorbed
    end do

    close(2)


end program pile

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
        ! Scatter angle
        call random_number(rn1)
        angle = acos((2. * rn1) - 1.)
    else
        ! Maximum recoil energy that can be imparted
        max_recoil_energy = ((4 * Z) / ((1 + Z)**2)) * energy_in
        call random_number(rn1)
        ! Probability distribution for recoil energies is flat from 0 to max
        max_recoil_energy = max_recoil_energy * rn1 
        energy_out = energy_in - max_recoil_energy
        ! Scatter angle
        eta = acos(sqrt((max_recoil_energy/energy_in)*(((1.+Z)**2)/(4.*Z))))
        angle = atan((sin(2.*eta))/((1./Z) - cos(2.*eta)))
    end if

    ! Update energy
    energy_in = energy_out 
end subroutine energy_loss


! init_energy: Samples the Watt distribution for Uranium-235 to choose a random
!              initial energy of a simulated neutron.
!   energy: the energy of the neutron (result in MeV)
subroutine init_energy(energy)
    implicit none 
    real, intent(out) :: energy 
    real :: K, L, M, x, y, rn1, rn2 
    real :: a, b 

    ! Values for Uranium-235
    a = 0.9
    b = 1.0

    K = 1 + (b / 8 * a)
    L = (K + sqrt(K**2 - 1)) / a 
    M = a*L - 1
    do while(.true.)
        call random_number(rn1)
        call random_number(rn2)
        x = -log(rn1)
        y = -log(rn2)
        if ((y - M*(x + 1))**2 < (b * L * x)) then 
            energy = L * x
            exit
        endif
    enddo
end subroutine init_energy

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
    real, intent(in) :: angle 
    real :: scatter_dir(3)

    real :: Ex, Ey, Ez

    real :: arg, rn1
    real :: alpha, beta, gamma, theta, phi
    real :: sco1, sco2
    real :: r(3,3) ! Rotation matrix
    real :: s0(3) 

    ! Normalize direction to a unit vector (in case it wasn't)
    in_dir = in_dir / norm2(in_dir)
    Ex = in_dir(1)
    Ey = in_dir(2)
    Ez = in_dir(3)

    beta = acos(Ez)

    ! Approximation
    if (abs(beta) < 0.027) then 
        alpha = 0.0
    else 
        arg = Ey / sin(beta)
        if (abs(arg) >= 1.0) then  ! huh??
            arg = arg / (1.0001 * abs(arg))
        end if 
        alpha = asin(arg)
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
    r(1,1) = cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma)
    r(1,2) = cos(beta)*sin(alpha)*cos(gamma)+cos(alpha)*sin(gamma)
    r(1,3) =-sin(beta)*cos(gamma)
    r(2,1) =-sin(gamma)*cos(beta)*cos(alpha)-sin(alpha)*cos(gamma)
    r(2,2) =-sin(gamma)*cos(beta)*sin(alpha)+cos(alpha)*cos(gamma)
    r(2,3) = sin(beta)*sin(gamma)
    r(3,1) = sin(beta)*cos(alpha)
    r(3,2) = sin(alpha)*sin(beta)
    r(3,3) = cos(beta)
    s0 = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    ! Unit propagation vector of the scattered particle in the original frame
    scatter_dir = matmul(r, s0)

    in_dir = scatter_dir

end subroutine euler

subroutine xsection_from_energy(energy, energies, xsections, xsect)
    real, intent(in) :: energy ! Input energy of incident neutron
    real, intent(in) :: energies(31882) ! Energy data
    real, intent(in) :: xsections(31882) ! Corresponding cross section data
    real, intent(out) :: xsect ! The resulting cross section for the energy

    integer :: index1, index2
    integer :: low, high, mid
    real :: x1, y1, x2, y2

    ! Initialize indices
    low = 1
    high = size(energies)

    ! Perform binary search to find nearest indices
    do while (low < high)
        mid = (low + high) / 2
        if (energy < energies(mid)) then
            high = mid 
        else 
            low = mid + 1
        end if
    end do 

    index1 = low - 1
    index2 = low

    x1 = energies(index1)
    y1 = xsections(index1)
    x2 = energies(index2)
    y2 = xsections(index2)

    ! Linear interpolation to find approx. cross section
    xsect = y1 + ((energy - x1) * (y2 - y1)) / (x2 - x1)

end subroutine xsection_from_energy

! Imposes toroidal boundary conditions to simulate infinite lattice
subroutine wrap_around(position, D)
    real, intent(inout) :: position(3)
    real, intent(in) :: D

    ! Wrap around x and y directions
    position(1) = modulo(position(1), 2.0 * D) - D
    position(2) = modulo(position(2), 2.0 * D) - D

    ! Wrap around z direction
    position(3) = modulo(position(3), 3.0 * D) - 1.5 * D
end subroutine wrap_around

subroutine run_simulation(num_neutrons, D,  num_absorbed, total_fission_xsects, &
                          energies, xsections)
    implicit none
    integer, intent(in) :: num_neutrons
    real, intent(in) :: D ! Graphite brick width (to be optimized)
    integer, intent(out) :: num_absorbed ! Number of absorbed neutrons
    real, intent(in) :: energies(31882) ! Energy values 
    real, intent(in) :: xsections(31882) ! Corresponding cross section data
    real, intent(out) :: total_fission_xsects ! Simulated total cross sections

    ! CONSTANTS 
    real :: NA   ! Avogadro's number
    real :: rho  ! Density of medium (graphite) 
    integer :: Z ! atomic number (carbon-12)
    real :: pi 
    real :: R    ! Uranium rod radius
    real :: R_squared ! rod radius squared

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
    real :: fission_xsect        ! Cross section of U-235 fission

    ! MISC
    integer :: i, j  ! For iterating
    real :: rn       ! For random numbers
    
    pi = 4 * atan(1.0)
    NA = 6.023e23 ! 1 mol
    rho = 2.16    ! g/cm^3
    Z = 12        ! g/mol
    R = 4.1275    ! cm (3.25 in diameter)
    R_squared = R**2 ! Just to calculate it only once

    total_fission_xsects = 0
    num_absorbed = 0

    ! Main simulation loop to do for each neutron
    simulation_loop: do i = 1, num_neutrons
        !print *,"Simulating neutron ",i
        call random_seed() 

        ! Initialize position to origin 
        position(:) = 0
        ! Reset fission cross section
        fission_xsect = 0

        ! Initialize energy 
        call init_energy(energy)

        ! Build unit velocity vector
        call random_number(phi)
        phi = 2*pi * phi
        call random_number(theta)
        theta = acos(2*theta - 1) ! Properly sampled 

       ! Convert angles to cartesian vector
        velocity = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]

        ! Individual neutron loop
        ! (Exits loop if neutron backscatters/traverses/is absorbed)
        neutron_loop: do while (.true.)

            ! If the neutron hasn't hit a U235 cylinder,
            ! Find cross sections for elastic scattering / absorption
            call cross_section(energy, xelastic, xabsorbed)

            ! Determine distance to next interaction
            call random_number(rn)
            ! Depends on total cross section
            dnext = -log(rn) * Z / ((xelastic + xabsorbed) * rho * NA)

            ! Calculate new position (at interaction)
            position = position + dnext * velocity

            ! Modular unit cell calculations
            ! "Wrap around" back into the unit cell
            call wrap_around(position, D)

            ! Check if neutron has struck a uranium rod
            ! Can be the rod it originated from or any other; doesn't matter
            !  since using toroidal geometry
            if ((position(1)**2 + position(2)**2 <= R_squared) .and. &
                (position(3) <= 1.5*D) .and. (position(3) >= -1.5*D)) then 
                call xsection_from_energy(energy, energies, xsections, fission_xsect)
                total_fission_xsects = total_fission_xsects + fission_xsect
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