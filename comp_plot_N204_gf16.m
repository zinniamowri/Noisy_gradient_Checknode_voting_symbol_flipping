
Eb_No_cnv=9:.5:11.5;
ber_cnv=[5.402e-03, 1.141e-03, 1.689e-04, 3.578e-05, 6.936e-06, 2.255e-06];

Eb_No_ngdsf=9:1:13;
uncoded= [1.804902e-02, 1.104902e-02, 6.125194e-03, 2.854495e-03, 1.046679e-03];
ber_ngdsf=[1.442157e-02,6.365196e-03, 1.749088e-03, 1.494500e-04,2.709264e-06];

%data from ems simulation
ems_Eb_N0=[8,8.4,9.2,9.6];
ems_ber=[0.0330501,0.00532898,0.000337068,4.8366e-06];


figure;

semilogy(Eb_No_ngdsf, uncoded, 'rx-', 'LineWidth', 1.2);
grid on;
hold on;

semilogy(Eb_No_cnv, ber_cnv, 'bx-','LineWidth', 1.2);

semilogy(Eb_No_ngdsf, ber_ngdsf, 'kx-','LineWidth', 1.2);

%semilogy(x_for_algo1, algorithm1, '-s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.2);


%semilogy(x_for_algo2, algorithm2, '-s', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2);

%semilogy(x_for_algo3, algorithm3, '-s','Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.2);

semilogy(ems_Eb_N0, ems_ber, 'mx-', 'LineWidth', 1.2);
hold off;

 ylim([10e-7 10e-2]);
 xlim([6 16]);
 %xticks([1 5 10 15 20]);
 xlabel('E_b/N_0 (dB)', 'FontSize', 14);
 ylabel('BER', 'FontSize', 14);
 %title('BER comparison for NB-LDPC code (45,2,3) over GF(8)');
 hold off;
 legend( '16-QAM (uncoded)', 'proposed NG-CNV-SF', 'NGDSF','EMS', 'Location', 'northeast', 'FontSize', 9);
