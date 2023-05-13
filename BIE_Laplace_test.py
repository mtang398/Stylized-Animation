import numpy as np
import scipy.linalg
import cmath
from time import time

tStart = 0.0
NBR_R = 10  # Resolution in radial direction
NBR_T = 10  # Resolution in boundary parameterisation
NBR_DOMAIN_POINTS = NBR_R * NBR_T  # Number of grid points in domain
NBR_PANELS = 30  # Number of GL-panels
NBR_POINTS_PER_PANEL = 16  # Number of points per GL-panel
NBR_PANEL_POINTS = NBR_POINTS_PER_PANEL * NBR_PANELS  # Useful abbrevation

gl_nodes = np.array([
    -0.989400934991650, -0.944575023073233, -0.865631202387832, -0.755404408355003,
    -0.617876244402644, -0.458016777657227, -0.281603550779259, -0.095012509837637,
    0.095012509837637, 0.281603550779259, 0.458016777657227, 0.617876244402644,
    0.755404408355003, 0.865631202387832, 0.944575023073233, 0.989400934991650
])

gl_weights = np.array([
    0.027152459411754, 0.062253523938648, 0.095158511682493, 0.124628971255534,
    0.149595988816577, 0.169156519395003, 0.182603415044924, 0.189450610455068,
    0.189450610455068, 0.182603415044924, 0.169156519395003, 0.149595988816577,
    0.124628971255534, 0.095158511682493, 0.062253523938648, 0.027152459411754
])

def create_grid(pz, NBR_R, NBR_T, tStart):
    for i in range(NBR_R):
        r = i * 0.999 / (NBR_R - 1)
        for j in range(NBR_T):
            t = 2.0 * np.pi * j / (NBR_T - 1)
            pz[i * NBR_T + j] = r * (1.0 + 0.3 * np.cos(5.0 * (t + tStart))) * np.exp(1j * (t + tStart))

def tau(ptau, t, N, tStart):
    for i in range(N):
        ptau[i] = (1.0 + 0.3 * np.cos(5.0 * (t[i] + tStart))) * np.exp(1j * (t[i] + tStart))

def taup(ptaup, t, N, tStart):
    for i in range(N):
        ptaup[i] = (-1.5 * np.sin(5.0 * (t[i] + tStart)) + 1j * (1.0 + 0.3 * np.cos(5 * (t[i] + tStart)))) * np.exp(1j * (t[i] + tStart))

def taupp(ptaupp, t, N, tStart):
    for i in range(N):
        ptaupp[i] = np.exp(1j * (t[i] + tStart)) * (-1.0 - 7.8 * np.cos(5.0 * (t[i] + tStart)) - (3.0 * 1j) * np.sin(5.0 * (t[i] + tStart)))

def gl16(pzDrops, pzDropsp, pzDropspp, pwDrops, ptpar, ppanels, NBR_PANELS, NBR_POINTS_PER_PANEL, gl_weights, gl_nodes):
    for i in range(NBR_PANELS):
        first = ppanels[i].real
        last = ppanels[i+1].real
        for j in range(16):
            pwDrops[i * NBR_POINTS_PER_PANEL + j] = 0.5 * (last - first) * gl_weights[j]
            ztmp = 0.5 * (first * (1 - gl_nodes[j]) + last * (1 + gl_nodes[j]))
            tau(pzDrops[i * NBR_POINTS_PER_PANEL + j:i * NBR_POINTS_PER_PANEL + j+1], [ztmp], 1, tStart)
            taup(pzDropsp[i * NBR_POINTS_PER_PANEL + j:i * NBR_POINTS_PER_PANEL + j+1], [ztmp], 1, tStart)
            taupp(pzDropspp[i * NBR_POINTS_PER_PANEL + j:i * NBR_POINTS_PER_PANEL + j+1], [ztmp], 1, tStart)

def init_domain(pz, pzDrops, pzDropsp, pzDropspp, ppanels, ptpar, pwDrops, NBR_R, NBR_T, NBR_PANELS, NBR_POINTS_PER_PANEL, gl_weights, gl_nodes, tStart):
    create_grid(pz, NBR_R, NBR_T, tStart)
    gl16(pzDrops, pzDropsp, pzDropspp, pwDrops, ptpar, ppanels, NBR_PANELS, NBR_POINTS_PER_PANEL, gl_weights, gl_nodes)
    for i in range(NBR_PANELS + 1):
        ppanels[i] = (1.0 + 0.3 * np.cos(5.0 * (ppanels[i] + tStart))) * np.exp(1j * (ppanels[i] + tStart))
    
def init_function(RHS, pu_ana, pzDrops, pz):
    zsrc1 = 1.5 + 1.5j
    zsrc2 = -0.25 + 1.5j
    zsrc3 = -0.5 - 1.5j

    for i in range(NBR_PANEL_POINTS):
        RHS[i] = np.real(1 / (pzDrops[i] - zsrc1) + 1 / (pzDrops[i] - zsrc2) + 1 / (pzDrops[i] - zsrc3))

    pumax = 0  # Used to obtain relative error.
    for i in range(NBR_DOMAIN_POINTS):
        pu_ana[i] = np.real(1 / (pz[i] - zsrc1) + 1 / (pz[i] - zsrc2) + 1 / (pz[i] - zsrc3))
        if abs(pumax) < abs(pu_ana[i]):
            pumax = abs(pu_ana[i])

    return pumax
    
def call_gsl_gmres(A_gmres, RHS, epsrel, compress):
    # Solve the linear system using scipy's linalg.solve function
    pmu = scipy.linalg.solve(A_gmres, RHS)
    return pmu

def solveDensity(pzDrops, pzDropsp, pzDropspp, pwDrops, RHS):
    NBR_PANEL_POINTS = len(pwDrops)

    A_gmres = np.zeros((NBR_PANEL_POINTS, NBR_PANEL_POINTS), dtype=np.complex128)

    for i in range(NBR_PANEL_POINTS):
        for j in range(NBR_PANEL_POINTS):
            A_gmres[i, j] = (1.0 / 2.0 + (1 / (2 * np.pi)) * pwDrops[i] * np.imag(pzDropspp[i] / (2 * pzDropsp[i]))) if i == j else ((1 / (2 * np.pi)) * pwDrops[j] * np.imag(pzDropsp[j] / (pzDrops[j] - pzDrops[i])))

    epsrel = 1.0e-14
    compress = 0

    pmu = call_gsl_gmres(A_gmres, RHS, epsrel, compress)

    return pmu
    
def computeSolution(pmu, pz, pwDrops, pzDrops, pzDropsp):
    NBR_R = int(np.sqrt(len(pz)))
    NBR_T = int(np.sqrt(len(pz)))
    NBR_PANEL_POINTS = len(pwDrops)
    
    pu = np.zeros((NBR_R, NBR_T))

    for i in range(NBR_R):
        for j in range(NBR_T):
            pu_sum = 0
            for k in range(NBR_PANEL_POINTS):
                pu_sum += pmu[k] * pwDrops[k] * np.imag(pzDropsp[k] / (pzDrops[k] - pz[i * NBR_T + j])) * 1.0 / (2.0 * np.pi)
            pu[i, j] = pu_sum

    return pu
    
def computeError(pu, pu_spec, pu_ana, pumax):
    NBR_DOMAIN_POINTS = len(pu)
    errorvec = np.zeros(NBR_DOMAIN_POINTS)
    errormax = 0

    for i in range(NBR_DOMAIN_POINTS):
        errorvec[i] = abs(pu_spec[i] - pu_ana[i]) / pumax
        if errormax < errorvec[i]:
            errormax = errorvec[i]

    print(f"Max error: {errormax:.16e}")

    return errorvec

# function for special laplace quadrature
startsize = 32

W16 = np.array([0.0271524594117541,
       0.0622535239386479,0.0951585116824928,0.1246289712555339,0.1495959888165767,
       0.1691565193950025,0.1826034150449236,0.1894506104550685,0.1894506104550685,
       0.1826034150449236,0.1691565193950025,0.1495959888165767,0.1246289712555339,
       0.0951585116824928,0.0622535239386479,0.0271524594117541])

W32 = np.array([0.0070186100094701,
       0.0162743947309057,0.0253920653092621,0.0342738629130214,0.0428358980222267,
       0.0509980592623762,0.0586840934785355,0.0658222227763618,0.0723457941088485,
       0.0781938957870703,0.0833119242269467,0.0876520930044038,0.0911738786957639,
       0.0938443990808046,0.0956387200792749,0.0965400885147278,0.0965400885147278,
       0.0956387200792749,0.0938443990808046,0.0911738786957639,0.0876520930044038,
       0.0833119242269467,0.0781938957870703,0.0723457941088485,0.0658222227763618,
       0.0586840934785355,0.0509980592623762,0.0428358980222267,0.0342738629130214,
               0.0253920653092621,0.0162743947309057,0.0070186100094701])

IP1 = np.array([0.7082336805923402,
			0.4174603456614055,0.1201306210679926,-0.0237384214490220,-0.0247115099296993,
			0.0087677652494516,0.0103484922682165,-0.0048829375742948,-0.0058171385997042,
			0.0033716100497611,0.0039063895362322,-0.0026884914506444,-0.0029826204379799,
			0.0023966133657277,0.0025278996374981,-0.0023520880120266,-0.3533989766039463,
			0.1289141029247849,0.4904381439593198,0.4476888554483156,0.1595815187430263,
			-0.0431996043883891,-0.0450542761365297,0.0198432780500969,0.0226533071187764,
			-0.0127661666173026,-0.0145091849486018,0.0098524148293410,0.0108283966165110,
			-0.0086459140673153,-0.0090837161649102,0.0084360245293210,0.2547888140892103,
			-0.0794386511305539,-0.1800216753917482,0.1127069946981099,0.4529903213227301,
			0.4571030934901456,0.1720456302558366,-0.0572686901163087,-0.0572070079705632,
			0.0298642753857761,0.0323480604597623,-0.0212795957958572,-0.0228944350026966,
			0.0180248076268620,0.0187756290920826,-0.0173659080015727,-0.1926427644762529,
			0.0575395983713891,0.1185502002165091,-0.0605587419142252,-0.1375504134801631,
			0.1105111980442395,0.4385063216111399,0.4620135345715323,0.1767505762588966,
			-0.0700000771781777,-0.0663627648029001,0.0404551610934199,0.0415343932539965,
			-0.0317585071715892,-0.0325161000341283,0.0298340633769956,0.1433033540000506,
			-0.0420553811514164,-0.0836233891757775,0.0401361142969071,0.0816535055732694,
			-0.0527563287484864,-0.1144994775853014,0.1112128815258087,0.4303680392263895,
			0.4662432580208928,0.1793331307963325,-0.0837301037551319,-0.0757157523283246,
			0.0540255806818779,0.0532626953367605,-0.0480491590014860,-0.0997192378404317,
			0.0290155272553322,0.0567456713515015,-0.0265072815309246,-0.0516841720948336,
			0.0312108723382504,0.0603950802150088,-0.0468713962117200,-0.0965418713375640,
			0.1113929985433533,0.4229333917674608,0.4726563404781661,0.1839456727063571,
			-0.1024390131967896,-0.0906825961622589,0.0783222643612123,0.0589566269168384,
			-0.0170808697411231,-0.0331317609494943,0.0152767692550160,0.0292158689453183,
			-0.0171477684792897,-0.0317987518605555,0.0230873922781534,0.0424630971897936,
			-0.0391477381829754,-0.0774309538666544,0.1069678364314479,0.4098502582609164,
			0.4887570918888832,0.2021622184713724,-0.1455583641763743,-0.0195214966778084,
			0.0056453278101818,0.0109121889216969,-0.0050042888041768,-0.0094951190796480,
			0.0055107724940780,0.0100569812321846,-0.0071340625232678,-0.0126690018860246,
			0.0110418399786724,0.0197819310583685,-0.0222335618307413,-0.0445659130687800,
			0.0796393408723433,0.3555539698235838,0.5967331669239306])

IP2 = np.array([0.7138621264850029,
			0.4158614649995460,0.1171390533885666,-0.0224309414525152,-0.0223867275212341,
			0.0075268332425387,0.0083097853751947,-0.0036134992927350,-0.0038983391485322,
			0.0020027759055358,0.0020013610561088,-0.0011449345392067,-0.0010004418238182,
			0.0005796227498290,0.0003691229777510,-0.0001148410895266,-0.3731117376310109,
			0.1345146978688941,0.5009196712113994,0.4431061809431234,0.1514292542409961,
			-0.0388453473755598,-0.0378952348465420,0.0153814075225089,0.0159014848426175,
			-0.0079431247886893,-0.0077862576814684,0.0043949156603830,0.0038044673660606,
			-0.0021902526754258,-0.0013893468075581,0.0004314370407234,0.2935334077534948,
			-0.0904491991507758,-0.2006375430165828,0.1217267282571176,0.4690505693866059,
			0.4485149826814497,0.1579049657964081,-0.0484399253991464,-0.0438186361102529,
			0.0202761928802727,0.0189425113789293,-0.0103579732562729,-0.0087773455062908,
			0.0049826169161143,0.0031336115858597,-0.0009691268935199,-0.2543216125483122,
			0.0750746089078695,0.1514059982920834,-0.0749489083470473,-0.1632097247818790,
			0.1242574593626827,0.4611915989585189,0.4478105302111243,0.1551400216478456,
			-0.0544610911859422,-0.0445314841467965,0.0225651764329818,0.0182471281387628,
			-0.0100600543577763,-0.0062187415151346,0.0019078707296151,0.2312943045160101,
			-0.0670850646237963,-0.1305709521800638,0.0607297941305177,0.1184505233059086,
			-0.0725218317815637,-0.1472268604148448,0.1317870430598982,0.4618288268307383,
			0.4434844549025628,0.1471232281426801,-0.0570984665414406,-0.0406678216286402,
			0.0209226990236611,0.0124538953296838,-0.0037566466272948,-0.2171239069846040,
			0.0624388430107100,0.1195285512638427,-0.0541067920843720,-0.1011439299315847,
			0.0578788932057221,0.1047623469874682,-0.0749282556026716,-0.1397580556688505,
			0.1429367300147700,0.4680721498192956,0.4348189017231927,0.1332828758133183,
			-0.0535184789189454,-0.0286039577327160,0.0082607580054626,0.2087875429056409,
			-0.0597829885222319,-0.1135080588890587,0.0507179130294734,0.0929917302083465,
			-0.0517207938386748,-0.0897133328585768,0.0600282764461917,0.0999806752076895,
			-0.0817026011454862,-0.1393794337995033,0.1600513710164597,0.4830068086060422,
			0.4153122181040704,0.1037159232533378,-0.0249698016066503,-0.2049002094488522,
			0.0585617629264828,0.1108029670564811,-0.0492413053477593,-0.0895742689267964,
			0.0492637410706939,0.0840953326985925,-0.0549762659931624,-0.0884105585949814,
			0.0683011464015245,0.1055383029995325,-0.0985990125544629,-0.1556639994547752,
			0.2005703021781758,0.5406401699812300,0.3033999036738149])

def vandernewton(T, b, N):
    for k in range(1, N):
        for j in range(N - 1, k - 1, -1):
            b[j] -= T[k - 1] * b[j - 1]

    for k in range(N - 1, 0, -1):
        for j in range(k, N):
            b[j] = b[j] / (T[j] - T[j - k])
            b[j - 1] -= b[j]

def vandernewtonT(T, b, N):
    for k in range(N - 1):
        for j in range(N - 1, k, -1):
            b[j] = (b[j] - b[j - 1]) / (T[j] - T[j - k - 1])

    for k in range(N - 1, -1, -1):
        for j in range(k, N - 1):
            b[j] -= T[k] * b[j + 1]

def IPmultR(in_, out):
    for i in range(16):
        t1 = 0
        t2 = 0
        ptr = i
        for j in range(8):
            t1 += IP1[ptr] * (in_[j] + in_[15 - j])
            t2 += IP2[ptr] * (in_[j] - in_[15 - j])
            ptr += 16

        out[i] = t1 + t2
        out[31-i] = t1 - t2

def specialquadlapl(u_specq, u_standardq, mu, zDom, zDrops, zpDrops, wDrops, panels):
    NBR_DOMAIN_POINTS = len(u_specq)
    NBR_PANELS = len(panels) - 1

    for i in range(NBR_DOMAIN_POINTS):
        u_specq[i] = u_standardq[i]
        zk = zDom[i]
        for k in range(NBR_PANELS):
            mid = (panels[k + 1] + panels[k]) / 2.0
            len_ = (panels[k + 1] - panels[k])

            if abs(zk - mid) < abs(len_):
                nz = 2 * (zk - mid) / len_
                oldsum = 0
                testsum = 0
                lg1 = cmath.log(1 - nz)
                lg2 = cmath.log(-1 - nz)

                tz = [zDrops[k * 16 + j] for j in range(16)]
                nzpan = [(2.0 * (tz[j] - mid) / len_) for j in range(16)]

                if -1 < nz.real < 1:
                    furthercheck = False
                    for j in range(16):
                        if (nz.imag > 0 and nzpan[j].imag > nz.imag) or (nz.imag < 0 and nzpan[j].imag < nz.imag):
                            furthercheck = True
                            break

                    if furthercheck:
                        tmpT = [nzpan[j].real for j in range(16)]
                        tmpb = [nzpan[j].imag for j in range(16)]

                        vandernewtonT(tmpT, tmpb, 16)
                        test = tmpb[15]

                        for j in range(14, -1, -1):
                            test = test * nz.real + tmpb[j]

                        if (nz.imag > 0 and test > nz.imag) or (nz.imag < 0 and test < nz.imag):
                            if nz.imag > 0:
                                lg1 -= np.pi * 1j
                                lg2 += np.pi * 1j
                            else:
                                lg1 += np.pi * 1j
                                lg2 -= np.pi * 1j

                p32 = np.zeros(33, dtype=complex)
                p32[0] = lg1 - lg2

                oldsum = 0
                testsum = 0
            for j in range(16):
                tzp = zpDrops[k * 16 + j]
                tf = mu[k * 16 + j]
                tW = wDrops[k * 16 + j]
                oldsum += tW * tf * (tzp.imag / (tz[j] - zk).imag) / (2. * np.pi)
                testsum += tW * tzp / (tz[j] - zk)

            if abs(p32[0] - testsum) > 1e-13:
                tf32 = np.zeros(32, dtype=complex)
                tz32 = np.zeros(32, dtype=complex)
                tzp32 = np.zeros(32, dtype=complex)

                IPmultR(tf, tf32)
                IPmultR(tz, tz32)
                IPmultR(tzp, tzp32)

                plen = tW[0] / W16[0]
                tW32 = np.zeros(32)
                orig32 = np.zeros(32, dtype=complex)
                o32sum = 0

                for j in range(32):
                    tW32[j] = W32[j] * plen
                    orig32[j] = tW32[j] / (tz32[j] - zk)
                    o32sum += tzp32[j] * orig32[j]

                if abs(o32sum - p32[0]) < 1e-13:
                    newsum = 0
                    for j in range(32):
                        newsum += tW32[j] * tf32[j] * (tzp32[j].imag / (tz32[j] - zk).imag) * 0.5 / np.pi

                    u_specq[i] += (newsum - oldsum).real

                else:
                    nzpan32 = np.zeros(32, dtype=complex)
                    IPmultR(nzpan, nzpan32)
                    sign = -1

                    for j in range(1, 33):
                        p32[j] = nz * p32[j - 1] + (1.0 - sign) / j
                        sign = -sign

                    vandernewton(nzpan32, p32, 32)

                    new1 = 0
                    for j in range(32):
                        new1 += (p32[j] * tf32[j]).imag * 0.5 / np.pi

                    modif = new1 - oldsum
                    u_specq[i] += modif.real

def main():
    umax = 0  # Used to obtain relative error.

    # Allocate memory for complex arrays
    pz = np.zeros(NBR_DOMAIN_POINTS, dtype=np.complex128)
    pzDrops = np.zeros(NBR_PANEL_POINTS, dtype=np.complex128)
    pzDropsp = np.zeros(NBR_PANEL_POINTS, dtype=np.complex128)
    pzDropspp = np.zeros(NBR_PANEL_POINTS, dtype=np.complex128)
    ppanels = np.zeros(NBR_PANELS + 1, dtype=np.complex128)

    # Allocate memory for double arrays
    RHS = np.zeros(NBR_PANEL_POINTS)
    ptpar = np.zeros(NBR_PANEL_POINTS)
    pwDrops = np.zeros(NBR_PANEL_POINTS)
    pmu = np.zeros(NBR_PANEL_POINTS)
    pu = np.zeros(NBR_DOMAIN_POINTS)
    pu_spec = np.zeros(NBR_DOMAIN_POINTS)
    pu_ana = np.zeros(NBR_DOMAIN_POINTS)
    perrorvec = np.zeros(NBR_DOMAIN_POINTS)

    time_start = time()

    # Initialize the domain
    init_domain(pz, pzDrops, pzDropsp, pzDropspp, ppanels, ptpar, pwDrops, NBR_R, NBR_T, NBR_PANELS, NBR_POINTS_PER_PANEL, gl_weights, gl_nodes, tStart)

    # Evaluate the given right-hand side and obtain the analytical solution
    init_function(RHS, pu_ana, pzDrops, pz, umax)

    time_init = time()

    # Solve for density pmu
    solveDensity(pzDrops, pzDropsp, pzDropspp, pwDrops, RHS, pmu)

    time_density = time()

    # Evaluate the solution pu
    computeSolution(pmu, pz, pwDrops, pzDrops, pzDropsp, pu)

    time_sol = time()

    # Evaluate the solution pu_spec with special quadrature
    specialquadlapl(pu_spec, pu, pmu, pz, pzDrops, pzDropsp, pwDrops, ppanels)

    time_spec_q = time()

    # Compute the error perrorvec
    computeError(perrorvec, pu, pu_spec, pu_ana, umax)

    print("Timings for run on starfish")
    print("Parameters:")
    print(f"Npanels = {NBR_PANELS}")
    print(f"NBR_R = {NBR_R} \t NBR_T = {NBR_T}")
    print(f"Ndomain_points = {NBR_DOMAIN_POINTS}")
    print(f"Time for density: {time_density - time_init}")
    print(f"Time for solution: {time_sol - time_density}")
    print(f"Time for special quad: {time_spec_q - time_sol}")
    print("\n")
    print(f"Total time: {time_spec_q - time_init}")

if __name__ == "__main__":
    main()