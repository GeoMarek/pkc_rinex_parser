package nav_proc;
import java.awt.*;
import java.io.*;
import java.util.StringTokenizer;


public class Stdalone {
    static final double We = 7.292115E-5;
    static final double c = 299792458.0;
    static final double pi = 3.1415926535898;
    static final double Wedot = 7.2921151467E-5;
    static final double mu = 3.986005E+14;
    static final double F = -4.442807633E-10;
    static final double a = 6378137.0;
    static final double b = 6356752.31;
    static final double e1sqr = (a * a - b * b) / (a * a);
    static final double e2sqr = (a * a - b * b) / (b * b);
    static final double ITERATIONS =9;

    /**
     * Convert LLA position to XYZ position and return it
     * @param Xi array of LLA coordinates
     * @return array of XYZ coordinates
     */
    public double[] LLA2XYZ(double[] Xi) {
        double N = a / Math.sqrt(1.0 - e1sqr * Math.sin(Xi[0]) * Math.sin(Xi[0]));
        double[] Xo = new double[3];
        Xo[0] = (N + Xi[2]) * Math.cos(Xi[0]) * Math.cos(Xi[1]);
        Xo[1] = (N + Xi[2]) * Math.cos(Xi[0]) * Math.sin(Xi[1]);
        Xo[2] = (N * (1.0 - e1sqr) + Xi[2]) * Math.sin(Xi[0]);
        System.out.println("jd LLA2XYZ = " +Xo[0] +" "+Xo[1]+" "+ Xo[2]+" R = " + Math.sqrt(Xo[0]*Xo[0]+Xo[1]*Xo[1]+Xo[2]*Xo[2]));
        return Xo;
    }

    public double[] XYZ2LLA(double[] Xi) {
        double p = Math.sqrt(Xi[0] * Xi[0] + Xi[1] * Xi[1]);
        double T = Math.atan((Xi[2] * a) / (p * b));
        double sT = Math.sin(T);
        double cT = Math.cos(T);
        double[] Xo = new double[3];
        Xo[0] = Math.atan((Xi[2] + e2sqr * b * sT * sT * sT) / (p - e1sqr * a * cT * cT * cT));
        if ( Xi[0] == 0.0 )
            Xo[1] = pi / 2.0 ;
        else
            Xo[1] = Math.atan(Xi[1] / Xi[0]);
        double N = a / Math.sqrt(1.0 - e1sqr * Math.sin(Xo[0]) * Math.sin(Xo[0]));
        Xo[2] = p / Math.cos(Xo[0]) - N;
        return Xo;
    }


    public double[] satpos(double[] eph, double Ttr) {

        double Crs = eph[0];
        double dn  = eph[1] ;
        double M0  = eph[2] ;
        double Cuc = eph[3];
        double ec  = eph[4];
        double Cus = eph[5];
        double A   = eph[6] * eph[6];
        double Toe = eph[7];
        double Cic = eph[8];
        double W0  = eph[9] ;
        double Cis = eph[10];
        double i0  = eph[11];
        double Crc = eph[12];
        double w   = eph[13];
        double Wdot= eph[14];
        double idot= eph[15];
        double T;
        T= Ttr - Toe;
        if ( T >  302400.0 ) T = T - 604800.0;
        if ( T < -302400.0 ) T = T + 604800.0;

        double n0 = Math.sqrt(mu / (A * A * A));
        double n = n0 + dn;

        double M = M0 + n * T;
        double E = M;
        System.out.println("jd  M mu MO T"+ " " +M + " " + mu + " " + M0 + " " + T);
        double Eold;
        do {
            Eold = E;
            E = M + ec * Math.sin(E);
        } while  ( Math.abs(E - Eold) >= 1.0e-8 );



        double snu = Math.sqrt(1 - ec * ec) * Math.sin(E) / (1 - ec * Math.cos(E));
        double cnu = (Math.cos(E) - ec) / (1 - ec * Math.cos(E));
        double nu;
        nu = Math.atan2(snu, cnu);

        double phi = nu + w;

        double du = Cuc * Math.cos(2 * phi) + Cus * Math.sin(2 * phi);
        double dr = Crc * Math.cos(2 * phi) + Crs * Math.sin(2 * phi);
        double di = Cic * Math.cos(2 * phi) + Cis * Math.sin(2 * phi);

        double u = phi + du;
        double r = A * (1 - ec * Math.cos(E)) + dr;
        double i = i0 + idot * T + di;

        double Xdash = r * Math.cos(u);
        double Ydash = r * Math.sin(u);

        double Wc= W0 + (Wdot - Wedot) * T - Wedot * Toe;

        double[] X = new double[3];
        X[0] = Xdash * Math.cos(Wc) - Ydash * Math.cos(i) * Math.sin(Wc);  //ECEF x
        X[1] = Xdash * Math.sin(Wc) + Ydash * Math.cos(i) * Math.cos(Wc);  //ECEF y
        X[2] = Ydash * Math.sin(i);  //ECEF z



        double Trel = F * ec * eph[6] * Math.sin(E); //rel

        return new double[] {X[0], X[1], X[2], Trel};
    }


    public double[] calcAzEl(double[] Xs, double[] Xu) {

        double x = Xu[0];
        double y = Xu[1];
        double z = Xu[2];
        double p = Math.sqrt(x * x + y * y);
        if ( p == 0.0 )
            return null;

        double R = Math.sqrt(x * x + y * y + z * z);
        //System.out.println("jd calcAzEl R = " + R);
        double[][] e = new double[3][3];
        e[0][0] = - y / p;
        e[0][1] = x / p;
        e[0][2] = 0.0;
        e[1][0] = - x * z / (p * R);
        e[1][1] = - y * z / (p * R);
        e[1][2] = p / R;
        e[2][0] = x / R;
        e[2][1] = y / R;
        e[2][2] = z / R;

        double[] d = new double[3];
        for (int k = 0; k < 3; k++) {
            d[k] = 0.0;
            for (int i = 0; i < 3; i++)
                d[k] += (Xs[i] - Xu[i]) * e[k][i];
        }

        double s = d[2] / Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);

        double El;
        if ( s == 1.0 )
            El = 0.5 * pi;
        else
            El = Math.atan(s / Math.sqrt(1.0 - s * s));

        double Az;
        if ( (d[1] == 0.0) && (d[0] > 0.0) )
            Az = 0.5 * pi;
        else if ((d[1] == 0.0) && (d[0] < 0.0) )
            Az = 1.5 * pi;
        else {
            Az = Math.atan(d[0] / d[1]);
            if ( d[1] < 0.0 )
                Az += pi;
            else if ( (d[1] > 0.0) && (d[0] < 0.0) )
                Az += 2.0 * pi;
        }

        return new double[] {Az, El};
    }


    public double ionocorr (double[] ion, double Latu, double Lonu, double Az, double El, double Ttr) {

        double a0 = ion[0];
        double a1 = ion[1];
        double a2 = ion[2];
        double a3 = ion[3];
        double b0 = ion[4];
        double b1 = ion[5];
        double b2 = ion[6];
        double b3 = ion[7];


        Latu = Latu / pi;
        Lonu = Lonu / pi;
        El = El / pi;


        double phi = 0.0137 / (El + 0.11) - 0.022;

        double Lati = Latu + phi * Math.cos (Az);
        if ( Lati > 0.416 )
            Lati = 0.416;
        else if ( Lati < -0.416 )
            Lati = -0.416;

        double Loni = Lonu + phi * Math.sin(Az) / Math.cos(Lati * pi);
        double Latm = Lati + 0.064 * Math.cos((Loni - 1.617) * pi);

        double T = 4.32E+4 * Loni + Ttr;
        if (T >= 86400.0 )
            T -= 86400.0;
        else if (T < 0 )
            T += 86400.0;

        double F = 1.0 + 16.0 * (0.53 - El) * (0.53 - El) * (0.53 - El);

        double per = b0 + b1 * Latm + b2 * Latm * Latm + b3 * Latm * Latm * Latm;
        if (per < 72000.0 )
            per = 72000.0;

        double x = 2 * pi * (T - 50400.0) / per;

        double amp = a0 + a1 * Latm + a2 * Latm * Latm + a3 * Latm * Latm * Latm;
        if (amp < 0.0 )
            amp = 0.0;

        double dTiono;
        if (Math.abs(x) >= 1.57 )
            dTiono = F * 5.0E-9;
        else
            dTiono = F * (5.0E-9 + amp * (1.0 - x * x / 2.0 + x * x * x * x / 24.0));

        return dTiono;
    }

    public double sub (double[][] A, int r, int c) {

        double[][] B = new double[3][3];
        int i1, j1;

        for (int i = 0; i < 3; i++) {
            i1 = i;
            if ( i >= r )
                i1++;
            for (int j = 0; j < 3; j++) {
                j1 = j;
                if ( j >= c )
                    j1++;
                B[i][j] = A[i1][j1];
            }
        }


        return B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
                - B[1][0] * (B[0][1] * B[2][2] - B[2][1] * B[0][2])
                + B[2][0] * (B[0][1] * B[1][2] - B[0][2] * B[1][1]);
    }


    public double[] solve(double[][] Xs, boolean[] SV,  double[] P, double[] Xr) {

        double[] R = new double[33];
        double[] L = new double[33];
        double [][] A = new double[33][4];
        double[] AL = new double[4];
        double [][] AA = new double[4][4];
        double [][] AAi = new double[4][4];
        double det;
        double[] D = new double[5];


        int it = 0;
        System.out.println("solve");
        do {
            it++;

            for (int prn = 1; prn <= 32; prn++)
                if ( SV[prn] ) {
                    System.out.println("jd1 X Y Z P "+ prn +" "+Xs[prn][0]+" "+ Xs[prn][1]+" "+ Xs[prn][2]+" "+P[prn]);

                    R[prn] =  Math.sqrt((Xr[0] - Xs[prn][0]) * (Xr[0] - Xs[prn][0])
                            + (Xr[1] - Xs[prn][1]) * (Xr[1] - Xs[prn][1])
                            + (Xr[2] - Xs[prn][2]) * (Xr[2] - Xs[prn][2]));

                    L[prn] = P[prn] - R[prn];

                    for (int k = 0; k < 3; k++)
                        A[prn][k] = (Xr[k] - Xs[prn][k]) / R[prn];
                    A[prn][3] = -1.0;
                }


            for (int k = 0; k <= 3; k++) {
                AL[k] = 0.0;
                for (int prn = 1; prn <= 32; prn++)
                    if ( SV[prn] )
                        AL[k] += A[prn][k] * L[prn];
            }


            for (int k = 0; k <= 3; k++)
                for (int i = 0; i <= 3; i++) {
                    AA[k][i] =0.0;
                    for (int prn = 1; prn <= 32; prn++)
                        if ( SV[prn] )
                            AA[k][i] += A[prn][k] * A[prn][i];
                }

            det = AA[0][0] * sub(AA,0,0) - AA[1][0] * sub(AA,1,0)
                    + AA[2][0] * sub(AA,2,0) - AA[3][0] * sub(AA,3,0);
            if ( det == 0.0 )
                return null;

            int j;
            int n;
            for (int k = 0; k <= 3; k++)
                for (int i = 0; i <= 3; i++) {
                    n = k + i;
                    if ( n % 2 != 0 )
                        j = -1;
                    else
                        j = 1;
                    AAi[k][i] = j * sub(AA,i,k) / det;
                }

            for (int k = 0; k <= 3; k++) {
                D[k] = 0.0;
                for (int i = 0; i <= 3; i++)
                    D[k] += AAi[k][i] * AL[i];
            }


            for (int k = 0; k < 3; k++)
                Xr[k] += D[k];

        } while ( (it < ITERATIONS)
                && ((Math.abs(D[0]) + Math.abs(D[1]) + Math.abs(D[2])) >= 1.0E-2) );

        double Cr = D[3];

        if ( it >= ITERATIONS ) {
            System.out.println("rozw it = " + it);
            return null;
        }

        return new double[] {Xr[0], Xr[1], Xr[2], Cr};
    }


    public void main() {

        double Trc = 0;
        double Cr;
        double[] Xlla = new double[3];
        double[] Xr;
        boolean[] SV = new boolean[33];
        double[][] Xs = new double[33][3];
        double[][] eph = new double[33][16];
        double[][] clk = new double[33][5];
        double[] ion = new double[8];
        double[] Praw = new double[33];
        double[] Pcor = new double[33];


        FileDialog d = new FileDialog(new Frame(), "Wybierz plik", FileDialog.LOAD);
        d.setFile("*");
        d.show();
        String file1 = d.getFile();
        if (file1 != null)
            file1 = d.getDirectory() + d.getFile();
        else return;

        BufferedReader br;
        try {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(file1)));
        }
        catch (Exception e) {
            System.out.println ("Exception = " + e);
            return;
        }

        try {
            String line;
            if ((line = br.readLine()) != null)
                System.out.println("komenda : " + line);
            if ((line = br.readLine()) != null)
                Trc = Double.parseDouble(line.trim());

            if ((line = br.readLine()) != null)
                System.out.println("komenda : " + line);
            for (int i = 0; i < 8; i++)
                if ((line = br.readLine()) != null)
                    ion[i] = Double.parseDouble(line.trim());

            if ((line = br.readLine()) != null)
                System.out.println("komenda : " + line);
            StringTokenizer st;
            for (int prn = 1; prn < 32; prn++)
                SV[prn] = false;
            int prn1 = 0;
            do {
                if ((line = br.readLine()) != null) {
                    st = new StringTokenizer(line, " \t");
                    prn1 = Integer.parseInt(st.nextToken());

                    if ( prn1 != 0 ) {
                        Praw[prn1] = Double.parseDouble(st.nextToken());
                        SV[prn1] = true;
                        System.out.println("jd prn1 Praw " + prn1 +" "+Praw[prn1]);
                    }
                }
                else line = br.readLine();
            } while (prn1 != 0);

            if ((line = br.readLine()) != null)
                System.out.println("komenda: " + line);
            while ((line = br.readLine()) != null) {
                prn1 = Integer.parseInt(line.trim());
                for (int i = 0; i < 16; i++) {
                    line = br.readLine();
                    eph[prn1][i] = Double.parseDouble(line.trim());
                    System.out.println("jd " + eph[prn1][i]);
                }
                System.out.println("End of efemerids");
                for (int i = 0; i <= 4; i++) {
                    line = br.readLine();
                    clk[prn1][i] = Double.parseDouble(line.trim());
                }
            }
            br.close();
        }
        catch (IOException e) {
            System.out.println ("IOException = " + e);
        }


        System.out.println("Start ");

        System.out.println("       54 18 0");
        Xlla[0] = 54.373810539 * pi / 180.0;
        Xlla[1] = 18.614489984 * pi / 180.0;
        Xlla[2] = 0;
        Xr = LLA2XYZ(Xlla);

        Cr = 0;
        for (int prn = 1; prn <= 32; prn++)
            Pcor[prn] = 0.075 * c;

        for (int pass = 1; pass <= 2; pass++) {
            System.out.println();
            System.out.println("-------------------------- PASS " + pass + " -------------------------");

            for (int prn = 1; prn <= 32; prn++)
                if ( SV[prn] ) {


                    double tau = (Pcor[prn] + Cr) / c;
                    double Ttr = Trc - tau;


                    double[] tmp4 = satpos(eph[prn], Ttr);

                    double alpha = tau * We;
                    Xs[prn][0] = tmp4[0] * Math.cos(alpha) + tmp4[1] * Math.sin(alpha);
                    Xs[prn][1] = - tmp4[0] * Math.sin(alpha) + tmp4[1] * Math.cos(alpha);
                    Xs[prn][2] = tmp4[2];
                    double Trel = tmp4[3];
                    System.out.println("jd SV     : " + prn + " "+ Xs[prn][0] + " " + Xs[prn][1] + " " + Xs[prn][2]);

                    double[] tmp3 = new double[3];
                    System.arraycopy(Xs[prn], 0, tmp3, 0, 3);
                    double[] tmp2 = calcAzEl(tmp3, Xr);
                    double Az, El;
                    if (tmp2 == null) {
                        System.out.println("Blad in calcAzEl - ");
                        return;
                    }
                    else {
                        Az = tmp2[0];
                        El = tmp2[1];
                    }
                    System.out.println("Az, El : " + prn + " " + (Az * 180.0 / pi) + " " + (El * 180.0 / pi) );

                    double dTclck = - clk[prn][0] + clk[prn][2] + clk[prn][3] * (Ttr - clk[prn][1])
                            + clk[prn][4] * (Ttr - clk[prn][1]) * (Ttr - clk[prn][1])
                            + Trel ;
                    System.out.println("dTclck= clk[0]) + clk[2]+ clk[3]*(Ttr - clk[1])+  clk[4] Trel:"
                            + " " +dTclck +"="+ (- clk[prn][0])+ " " +clk[prn][2]+ " " + clk[prn][3]+ " " + (Ttr - clk[prn][1]) + " "
                            + clk[prn][4] + " " +  Trel);

                    double dTiono = ionocorr(ion, Xlla[0], Xlla[1], Az, El, Ttr);

                    double dRtrop = 2.312 / Math.sin(Math.sqrt(El * El + 1.904E-3))
                            + 0.084 / Math.sin(Math.sqrt(El * El + 0.6854E-3));

                    System.out.println("Corr   : " + prn + " " + Praw[prn]+ " " + (dTclck * c) + " " + (dTiono * c) + " " + dRtrop);

                    Pcor[prn] = Praw[prn] + dTclck * c - dTiono * c - dRtrop + Cr;

                }

            double[] tmp4 = solve(Xs, SV, Pcor, Xr);
            if (tmp4 == null) {
                System.out.println("Esprawdz");
                return;
            }
            else {
                Xr[0] = tmp4[0];
                Xr[1] = tmp4[1];
                Xr[2] = tmp4[2];
                Cr = tmp4[3];
            }
            System.out.println();
            System.out.println("Pos XYZ: " + Xr[0] + " " + Xr[1] + " " + Xr[2] + " " + Cr);
            Xlla = XYZ2LLA(Xr);
            System.out.println("Pos LLA: " + (Xlla[0] * 180.0 / pi) + " " + (Xlla[1] * 180.0 / pi) + " " + Xlla[2]);
        }

    }
}
