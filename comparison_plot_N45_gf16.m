%BER Performance of NG-CNV-SF over GF(16) with 16-QAM, N=45, K=15
%nbldpc_45_2_3_gf16_seed7.mat
%algorithm1 (Paper name: )
%EMS
%simulation setup of our own

Eb_No_db=[9.5, 10.5, 11.5, 12.5, 13.5, 14.5];

uncoded=[6.313838e-02, 4.460620e-02, 2.700759e-02, 1.526958e-02, 7.541333e-03, 3.197833e-03];

CNV=[2.063416e-02, 6.979734e-3, 1.628819e-3, 3.375400e-4, 3.766667e-5, 3.333333e-6];

NGDSF = [3.507218e-02, 1.457211e-02, 4.780164e-03, 2.125148e-03, 5.245214e-04, 1.547988e-04];

%algorithm1=[1.709049e-02, 5.612815e-03, 1.827029e-03, 4.800646e-04, 1.035000e-04,2.033333e-05];

%data from plotdizitizer
x_for_algo1= [10.2, 11, 12, 12.5, 13.5, 14.5, 15.5];
algorithm1=[0.0040, 0.0016, 0.00036, 0.000072, 0.00001199, 0.00000188, 1.35e-7];

%algorithm2=[ 0.01167, 0.00344, 0.000951, 0.0002241, 1.1e-05, 6.333333e-06]

%data from plotdizitizer
x_for_algo2=[6.2, 6.8, 7, 7.4, 7.8, 8, 8.4, 8.9, 9.4, 9.5];
algorithm2=[0.00103565, 0.0005800, 0.000119, 0.000069776, 0.000034189, 0.0000143315, 0.00000269, 0.00000129, 3.3e-7, 1.39e-8];

x_for_algo3=[9.6, 10.2, 11.17, 11.96, 12.7, 13.5];
algorithm3=[0.0003473, 0.0001089793, 0.0000273569, 0.00000439697, 0.000001032360, 1.1614e-7];
Eb_No_ems=[3.6, 3.8, 4, 4.2, 4.4];
ber_ems=[3.4327e-05, 2.12528e-05, 7.33333e-06, 4.66667e-06, 1.33333e-06];


figure;

semilogy(Eb_No_db, uncoded, 'rx-', 'LineWidth', 1.2);
grid on;
hold on;

semilogy(Eb_No_db, CNV, 'bx-','LineWidth', 1.2);

semilogy(Eb_No_db, NGDSF, 'kx-','LineWidth', 1.2);

semilogy(x_for_algo1, algorithm1, '-s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.2);


semilogy(x_for_algo2, algorithm2, '-s', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2);

semilogy(x_for_algo3, algorithm3, '-s','Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.2);

semilogy(Eb_No_ems, ber_ems, 'mx-', 'LineWidth', 1.2);
hold off;

 ylim([10e-9 10e-1]);
 xlim([1 20]);
 xticks([1 5 10 15 20]);
 xlabel('E_b/N_0 (dB)', 'FontSize', 14);
 ylabel('BER', 'FontSize', 14);
 %title('BER comparison for NB-LDPC code (45,2,3) over GF(16)');
 hold off;
 legend( 'uncoded (normalized Eb/No)', 'proposed NG-CNV-SF', 'NGDSF', 'algorithm1','algorithm2','algorithm3', 'EMS', 'Location', 'northeast', 'fontsize', 6);
