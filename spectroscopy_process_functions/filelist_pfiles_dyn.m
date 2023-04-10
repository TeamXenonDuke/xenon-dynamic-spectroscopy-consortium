% All spectroscopy Pfiles

% IPF = [1:17];
% healthy = [18:29,35:46];
PAH = [54:55];

IPF = [1,5,7,9,11,13,15,16,17]; % stats IPF = [2,5,7,9,11,13,15,16,17]; % stats
healthy = [21,22,27,28,29,35,41,42]; %stats [19,20,23,25,27,28,29,35,41,42];


names = {'46B','46C','67','67A','67B','68','68A','69','69B','70','73A',...
    '76','76A','79','79A','80','003-007','3B','6A','39B','41A','42A','42B',...
    '49A','49D','55','64','65','003-001','93s1','93s2','93s3','93s4',...
    '93s5','93A','94s1','94s2','94s3','94s4','94s5','94A','95s1','95s2',...
    '95s3','95s4','95s5','59','74','87','88','91','66','66A'}';

FVC = [2.0600    1.7800    2.5400    2.5100    2.5100    1.9300    1.8600    4.4300    4.3000    3.100,...
    2.6600    1.6800    1.6100    2.6700    3.0400    3.5500    2.3100    4.5200    5.3600    4.3900,...
    5.6200    6.1600    6.1600    4.1700    3.7400    2.7700    5.0300    3.1800    3.2900         0,...
    0     0     0     0     0     0     0     0     0     0 6.3300         0         0         0,...
    0         0    4.9100    2.4600         0         0     0     0     0];


DLco = [7.0000    6.1000    9.7000   11.9000   11.6000   11.5000   11.6000   18.0000   17.0000,...
    12.8000   12.8000    8.9000    8.9000   14.8000   14.6000   13.6000   11.4000         0          0,...
    25.1000         0         0         0         0   27.9000         0   28.8000   26.5000   16.8000,...
    0         0         0         0         0         0         0         0         0         0         0,...
    0         0         0         0         0         0         0         0         0         0         0         0         0];
    
loc = 'D:\OneDrive\Documents\Duke\CIVM Research\pfiles\';

all_pfiles_dyn = {

    'Subject_002_046B\P17920.7' % IPF
    'Subject_002_046\ScanC\P_002-046C_BH_SPECTROSCOPY_19456.7' % IPF
    'Subject_002_067\P06144.7' % IPF
    'Subject_002_067A\P09728.7' % IPF
    'Subject_002_067B\P15872.7' % IPF    
    'Subject_002_068\P16384.7' % IPF
    'Subject_002_068A\P_002-068A_BH_SPECTROSCOPY_14848.7' % IPF
    'Subject_002_069\P19456.7' % IPF    
    'Subject_002_069B\P18944.7' % IPF new
    'Subject_002_070\P_002-070_BH_SPECTROSCOPY_12800.7' % IPF
    'Subject_002_073A\P12288.7' % IPF new
    'Subject_002_076\P07168.7' % IPF
    'Subject_002_076A\P43008.7' % IPF
    'Subject_002_079\P_002-079_BH_SPECTROSCOPY_14848.7' % IPF
    'Subject_002_079A\P13312.7' % IPF new
    'Subject_002_080\P_002-080_BH_SPECTROSCOPY_17408.7' % IPF
    'Subject_003_007\P10240.7' % IPF new
    
    'Subject_002_003\ScanB\Series7_SPECTROSCOPY\P05120.7' % Healthy
    'Subject_002_006\ScanA\Series7_SPECTRO_2_PY_1\P49664.7' % Healthy
    'Subject_002_039\ScanB\Series4_CALIBRATION\P08192.7' % Healthy new
    'Subject_002_041\ScanA\Series7_SPEC_300_OPY\P11776.7' % Healthy
    'Subject_002_042\ScanA\Series7_SPECT_DYN_PY\P03072.7' % Healthy new
    'Subject_002_042\ScanB\Series7_SPEC_DYN_OPY\P04608.7' % Healthy new
    'Subject_002_049\ScanA\Dynamic_Spec\P09216.7' % Healthy new
    'Subject_002_049\ScanD\Series7_CALIBRATION\P17408.7' % Healthy
    'Subject_002_055\Series3_CALIBRATION\P14336.7' % Healthy
    'Subject_002_064\P03584.7' % Healthy
    'Subject_002_065\P35840.7' % Healthy
    'Subject_003_001\P09728.7' % Healthy
   
    'Subject_002_093\s1\P19968.7'
    'Subject_002_093\s2\P20480.7'
    'Subject_002_093\s3\P20992.7'
    'Subject_002_093\s4\P21504.7'
    'Subject_002_093\s5\P22016.7'
    'Subject_002_093A\P25600.7'
   
    'Subject_002_094\s1\P02560.7'
    'Subject_002_094\s2\P03584.7'
    'Subject_002_094\s3\P04608.7'
    'Subject_002_094\s4\P05632.7'
    'Subject_002_094\s5\P06656.7'    
    'Subject_002_094A\P22016.7'

    'Subject_002_095\s1\P16384.7'
    'Subject_002_095\s2\P16896.7'
    'Subject_002_095\s3\P17408.7'
    'Subject_002_095\s4\P17920.7'
    'Subject_002_095\s5\P18432.7'
    
    'Subject_002_059\Series8_CALIBRATION\P07680.7' % PAH
    'Subject_002_074\P06656.7' % PAH
    
    'Subject_002_087\P15360.7' % Emphy
    'Subject_002_088\P09216.7' % Mystery
    'Subject_002_091\P10752.7' % Mystery
    'Subject_002_066\P16384.7' % Pre Stent
    'Subject_002_066A\P14336.7' % Post Stent  

    }; 

all_pfiles_dyn = strcat(loc,all_pfiles_dyn);

%% Un-used pfiles below


%     'D:\Elly\Documents\Duke\CIVM Research\pfiles\Subject_002_069A\P14848.7' % Corrupt, pfile only has 32 frames w/ 2048 pts
%     'D:\Elly\Documents\Duke\CIVM Research\pfiles\Subject_002_041\ScanB\Series7_SPEC_DYN_OPY\P19456.7' % Healthy new

% SUBJECT 92 has too low of SNR due to underwire bra
%     'D:\Elly\Documents\Duke\CIVM Research\pfiles\Subject_002_092\s1\P10752.7' % Healthy 
%     'D:\Elly\Documents\Duke\CIVM Research\pfiles\Subject_002_092\s2\P11264.7' % Healthy
%     'D:\Elly\Documents\Duke\CIVM Research\pfiles\Subject_002_092\s3\P12288.7' % Healthy long
%     'D:\Elly\Documents\Duke\CIVM Research\pfiles\Subject_002_092\s4\P12800.7' % Healhty long
%     'D:\Elly\Documents\Duke\CIVM Research\pfiles\Subject_002_092\s5\P13312.7' % Healhty exercise

