function Fixation_Resize(sample_rate_GT, sample_rate_ASL)

%Make sure arguments exist, if not add values (Weird non-transfer of
%arguments from one function to the other)
if nargin == 1 || isempty(sample_rate_GT) % If no sample rate for GT was input, this is default value
    sample_rate_GT = (1/60);
    disp('Missing argument');
elseif nargin < 1 || isempty(sample_rate_ASL) % If we are missing both arguments
    sample_rate_ASL =(1/30);
    sample_rate_GT = (1/60);
    disp('No valid Sample-Rate ASL Input');
elseif nargin == 2 %Everything's fine
    disp('GT sample rate is');
    disp(sample_rate_GT);
    disp('ASL sample rate is');
    disp(sample_rate_ASL);
    disp('Rawwwwr');
end


%Headers for each Page
header_1 = {'Left AOI', 'Start Time', 'X Position', 'Y Position', 'Frame Count', 'Stop_Time', 'GT Time', 'Compare',   'AOI #', 'AOI Name', '#of Fixations', 'Average Fixation Duration', '', }; %First AOI column header
header_2 = {'Middle AOI', 'Start Time', 'X Position', 'Y Position', 'Frame Count','Stop Time', 'GT Time', 'Compare',  'AOI #', 'AOI Name', }; % Second AOI column header
header_3 = {'Right AOI',  'Start Time', 'X Position', 'Y Position', 'Frame Count', 'Stop Time', 'GT Time', 'Compare', 'AOI #', 'AOI Name', }; % Third AOI column header

%ReadIn for the data
Start_Time = xlsread('Take_4_Oct_18' , 'J:J');  %Time the fixation begun
horz_posse = xlsread('Take_4_Oct_18' , 'Q:Q');  %x position
vert_posse = xlsread('Take_4_Oct_18' , 'R:R');  %y position
fix_frames = xlsread('Take_4_Oct_18' , 'K:K');  %The number of times the value is repeated (duration of fixation)
Stop_Time = xlsread('Take_4_Oct_18' , 'M:M');   %Time when the fixation stopped
Start_Frame = xlsread('Take_4_Oct_18' , 'T:T'); %Start Frame
Stop_Frame = xlsread('Take_4_Oct_18' , 'U:U');  %Stop Frame
Inter_Fix = xlsread('Take_4_Oct_18' , 'N:N');   %This is for making sure that the time between frames is correct (not a negative value)
AOI_Desg = xlsread('Take_4_Oct_18' , 'H:H');    %Based on the # of the AOI (Outside = 0, other AOIs are from 1 to i (i being the # of AOIs)
[num,txt,bling] = xlsread('Take_4_Oct_18', 'I:I');     %Baed on the name associated with the AOI (left, middle, and center)

[row_num, Col_num] = size(vert_posse);
generate = row_num; %generate now has the original number of rows in the document (interchanagable with row_num)

%Loop responsible for loading all text values (AOI names)
for xx=1:generate
    if xx == 1
      disp('Initializing AOI Names'); 
    
    else
    AOI_Name{xx} = (txt(xx));
    %disp(AOI_Name{xx});
    %disp('?');
    end
    
end
vertcat(AOI_Name);

%% Here goes function for catching Interfix duration problems 
% Basic idea is that it goes through all of the values in the
%'InterfixDuration' column and makes sure that they're all >0 (therefore
% possible)
for g=1:generate
    if Inter_Fix(g) < 0
        
        disp('White Whale detected Captain!');
        
        %Now removing values associated with defunct data 
        Start_Time(g,:) = [];
        horz_posse(g,:) = [];
        vert_posse(g,:) = [];
        fix_frames(g,:) = [];
        Stop_Time (g,:) = []; %
        Start_Frame(g,:)= []; %
        Stop_Frame(g,:) = [];
        AOI_Desg(g,:)   = [];
        AOI_Name{g,:}   = [];
    else
        %disp('Clear');
        %disp(AOI_Desg(g));
    end
    
end
%%

%Recalculating 'generate' (because it is used later on)
[row_num, Col_num] = size(vert_posse);
generate = row_num;

%Alters the values of fix_frames to be 60/second as opposed to 30/second.
for u=1:generate
    fix_frames(u) = (fix_frames(u)/sample_rate_GT); % duration multiplied by Sample rate of GT
end

%This is where the avergae fixation duration is calculated
steve=0;
for u=1:generate
    
    steve_1 = fix_frames(u,1);
 %   disp(steve_1);
    steve = steve + steve_1;  % Add all number together... to create ULTRA-STEVE
 %   disp(steve);
end

%Final Averaged Duration of the Fixations
Avg_Duartion = (steve/generate) * sample_rate_GT;

%%
%X & Y values for first AOI (as well as other pertinent info)
vert_posse_1 = []; % Y value
hori_posse_1 = []; % X value
Start_time_1 = []; % Start Time
fix_frames_1 = []; % Number of frames in the fixation
St1_time     = []; % Stop Time (Used with Start time to determine how many frames take place between fixations)
AOI_Identity = []; % ASLR's assigned AOI # associated with fixation
AOI_Namer    = {}; % String associated with the AOI
%

%X & Y values for 2nd AOI
vert_posse_2 = [];
hori_posse_2 = []; 
%

%X & Y values for 3rd AOI
vert_posse_3 = [];
hori_posse_3 = [];
%

%Important variables used in the 
TIME = [];                          % For comparason loop
AOI_NUM = [];                       % Identifies the AOI of the fixation
Compare = [];                       % This is for comparing TIME to 'Start time'
Holder=0;                           % Iterator value for debugging
Offical_S = (Start_Time(1,1)-1/60); % Starting value 
c=0;                                % Iterator value for debugging
Stringy = 0;                        % For checking the AOI #
%

%% Main for loop in which the ASLR data is converted to GT format
disp('entering for_loop');

%Loop through n times, where n is the # of rows in the excel document
for i = 1:generate
    
    %Take the number of frames the fixation lasts for.
    
    %x matrix of length 'fixation frame #'
    Calc_x_row = [fix_frames(i), 1];
    %y matrix of length 'fixation frame #'
    Calc_y_row = [fix_frames(i), 1];
    %Start_time Matrix
    Calc_S_row = [fix_frames(i), 1];
    %frame_numb Matrix
    Calc_f_row = [fix_frames(i), 1];
    %now we add the values into the rows
    Calc_St_row = [fix_frames(i), 1];
    %AOI # matrix
    Calc_Mnum_row =[fix_frames(i),1];
    %AOI Name
    Calc_Name_row ={fix_frames(i),1};
    %disp(Calc_Name_row);
    
  for x = 1:fix_frames(i) %extends the x/y/start values as many times as their are frames 

    if((horz_posse(i) > 1) || (vert_posse(i) > 1))  %If X or Y is larger than 1, then we know that it's not in the AOI (Fixation does not occur in AOI, and therefore must of occured in Outside)
        
    horz_posse(i)=  29; % Default value meaning 'GT is going to do nothing with this data
    vert_posse(i)= -54; % Second verse, same as the first
    
    %Iterate the GT time
    Offical_S = Offical_S + sample_rate_GT;
    
    Calc_x_row(x,1) = horz_posse(i);                                % for all rows make the value x 
    Calc_y_row(x,1) = vert_posse(i);                                % for all rows make the value y
    Calc_S_row(x,1) = Start_Time(i)+(sample_rate_GT*(x-1));         % for all rows make the value start time + add sample_rate
    Calc_f_row(x,1) = fix_frames(i);                                % for all rows make the value frame #
    Calc_St_row(x,1) = Stop_Time(i);                                % End Time
    Calc_Mnum_row(x,1) = AOI_Desg(i);                               % Which AOI it fell in (Always 0 in this case)
    Calc_Name_row{x,1} = 'Outside';
    
    Holder = Calc_S_row(x,1); % Give holder most recent value of start time used as a comparason value for removing values
    
    
    else  %This means the fixation fell in a lookzone   

    %Iterate the GT time
    Offical_S = Offical_S + sample_rate_GT;
    

        %Here, regardless of the AOI fixated in, the necessary
        %manipulations are applied.  The output for each monitor can be
        %altered below.


        horzi_posse(i) = horz_posse(i) * 1440; % alter the value to match the x value (horz_posse)
        verti_posse(i) = vert_posse(i) * 900;  % alter the value to match the y value (vert_posse)   
        
        Calc_x_row(x,1) = horzi_posse(i);                           % for all rows make the value x 
        Calc_y_row(x,1) = verti_posse(i);                           % for all rows make the value y
        Calc_S_row(x,1) = Start_Time(i)+(sample_rate_GT*(x-1));     % for all rows make the value start time + add sample_rate
        Calc_f_row(x,1) = fix_frames(i);                            % for all rows make the value frame #
        Calc_St_row(x,1) = Stop_Time(i);                            % End Time
        Calc_Mnum_row(x,1) = AOI_Desg(i);                           % Which AOI it fell into
        
        
  %This series of if statements is present to assign the 'name' associated with the AOI 
  %Unfortunately I couldn't find a more simple way to accomplish this
    if     strcmpi('Left AOI', AOI_Name{i})
               Calc_Name_row{x,1} = 'Left AOI';
                
    elseif strcmpi('Middle AOI', AOI_Name{i})
               Calc_Name_row{x,1} = 'Middle AOI';
            
    elseif strcmpi('Right AOI', AOI_Name{i})
               Calc_Name_row{x,1} = 'Right AOI';
                
    else
               Calc_Name_row{x,1} = AOI_Name{i};

    end
     
   
   Holder = Calc_S_row(x,1); %Start value for checking how far apart projected time is from ASLR time
   
    end
    
  
  end
        
       % One last concatination
       hori_posse_1 = vertcat(hori_posse_1, Calc_x_row); 
       vert_posse_1 = vertcat(vert_posse_1, Calc_y_row);
       Start_time_1 = vertcat(Start_time_1, Calc_S_row);
       fix_frames_1 = vertcat(fix_frames_1, Calc_f_row);
       St1_time     = vertcat(St1_time, Calc_St_row); 
       AOI_Identity = vertcat(AOI_Identity, Calc_Mnum_row);
       AOI_Namer    = vertcat(AOI_Namer, Calc_Name_row);
       
       
    %%Okay here I'm going for adding space in between the fixations 
    %It should be equal to the Start time of i+1 - stop time of i
    
    if(i < row_num)
    the_dude = Start_Time(i+1);
    the_anti = Stop_Time(i);    
    
    extra_space =(the_dude-(the_anti));
    extra_space =((extra_space/(sample_rate_GT)));
    %Okay so I can do a time check here, where I ask the program the time
    %I've displayed is equilvalent to the start time of the fixation

    
    %Arrays for the concatinating 
    Ex_x = [extra_space, 1];
    Ex_y = [extra_space, 1];
    Ex_S = [extra_space, 1];
    Ex_f = [extra_space, 1];
    Ex_St = [extra_space, 1];
    Ex_AN = [extra_space, 1];
    Ex_NM = {extra_space, 1};

        %Add the values to the array in this loop
        for y=1:(extra_space+1)
        
        %Iterator values being updated
        Offical_S = Offical_S + sample_rate_GT;
        Holder = Holder + sample_rate_GT;
        
        Ex_y(y,1) = -54;
        Ex_x(y,1) = 29;
        Ex_S(y,1) = 0;
        Ex_f(y,1) = 0;
        Ex_St(y,1) = 0;
        Ex_AN(y,1) = 0;
        Ex_NM{y,1} = 'None';
        
        end
       
       %Vertical concatination for the first AOI 
       hori_posse_1 = vertcat(hori_posse_1, Ex_x);
       vert_posse_1 = vertcat(vert_posse_1, Ex_y);
       Start_time_1 = vertcat(Start_time_1, Ex_S);
       fix_frames_1 = vertcat(fix_frames_1, Ex_f);
       St1_time     = vertcat(St1_time, Ex_St); 
       AOI_Identity = vertcat(AOI_Identity, Ex_AN);
       AOI_Namer    = vertcat(AOI_Namer, Ex_NM);
          
    end
    
    %Is the projected GT time larger (by more than 1 frame) than the time
    %according to ASLR+
    %If so, delete the last row
    if  (Offical_S - Holder) >= (sample_rate_GT)

    disp('Time has exceeded Bounds set by ASLR');
    Offical_S = Offical_S - sample_rate_GT; %Decrease the GT time by one frame
    disp('performing incision');
    
    
    c = c + 1; % how many times has this operation been performed
    disp(c);
    
    %Get the most recent frame #
    [choka, choka_2] = size(fix_frames_1);
    
    hori_posse_1(choka,:) = []; % delete last x value
    vert_posse_1(choka,:) = []; % delete last y value
    Start_time_1(choka,:) = []; % Delete last Start time value
    fix_frames_1(choka,:) = []; % Deletes last Frame value
    St1_time(choka,:)     = []; % Deletes last St1 value
    AOI_Identity(choka,:) = []; % Deletes last AOI value
    AOI_Namer(choka,:)     = []; % Deletes last AOI Name
    

    %This is for the case that GT time has fallen behind the time according
    %to ASLR, why, I have no idea.
    
    elseif (Offical_S + sample_rate_GT) < (Holder)
    [choka, choka_2] = size(fix_frames_1);
    disp ('Time is lagging behind ASLR Data');
    disp (Offical_S);
    disp (Holder);
    
    end
    
end
%%
%For the check and the the intialization of the the time (first second is
%the time of the initial fixation)
[useful,not] = size(St1_time);

ac_time = Start_time_1(1,1); % ac_time now has the initial value of the first fixation
disp(useful);                % How many total rows there are

%For loop dealing with the comparason time
for T = 1:useful
    
     TIME(T,1) = ac_time;
     
     % The whole point of this is to compare generated time with time
     % according to the ASLR+
     check_1 = Start_time_1(T,1);
     
     if (check_1 == 0) % if the check takes place during space inbetween fixations
     Compare(T,1) = 0;
         
     else    % if the check takes place where there is time dats
     Transfer = check_1- ac_time;
     Compare(T,1) = Transfer;
     
     end
     ac_time = (ac_time + (sample_rate_GT)); % Updates the time through each loop
     
end

%%
%Here, I copy first AOI array, and then sort the fixation values based on
%which AOI they occur in

%disp(AOI_Namer);

[Array_2,useless] = size(AOI_Namer);
disp(Array_2);

hori_posse_2 = hori_posse_1; % Give X array to AOI 2
vert_posse_2 = vert_posse_1; % Give Y array to AOI 2
hori_posse_3 = hori_posse_1; % Give X array to AOI 3
vert_posse_3 = vert_posse_1; % Givy Y array to AOI 3

% Okay, basic concept here is pretty simple/capable of being expanded upon 
% Copy the two x/y arrays, in this loop we go through the array and change
% values based on which AOI the fixations occured in
for tt = 1:Array_2
    
   
    %Assign to 1st Matrix (nullify 2nd and 3rd AOIs)
    if strcmpi(AOI_Namer{tt}, 'Left AOI')
    %   disp('Left AOI');
        vert_posse_2(tt) = -54; %Y position null second column
        hori_posse_2(tt) =  29; %X position null second column
        vert_posse_3(tt) = -54; %Y position null third column
        hori_posse_3(tt) =  29; %X position null third column
     
    %Assign to 2nd Matrix (nullify 1st and 3rd AOIs)    
    elseif strcmpi(AOI_Namer{tt}, 'Middle AOI') 
    %   disp('Middle AOI');
        vert_posse_1(tt) = -54; %Make the first column Y null
        hori_posse_1(tt) =  29; %Make the first column X null
        vert_posse_3(tt) = -54; %Make the third column Y null
        hori_posse_3(tt) =  29; %Make the third column X null
        
    %Assign to 3rd Matrix (nullify 1st and 2nd AOIs)    
    elseif strcmpi(AOI_Namer{tt}, 'Right AOI')
    %   disp('Right AOI');
        vert_posse_1(tt) = -54; %Make the first column Y null
        hori_posse_1(tt) =  29; %Make the first column X null
        vert_posse_2(tt) = -54; %Make the second column Y null
        hori_posse_2(tt) =  29; %Make the second column X null
        
    end
    
   
    vertcat(hori_posse_2);
    vertcat(vert_posse_2);
    
    vertcat(hori_posse_3);
    vertcat(vert_posse_3);
    
end
%%
disp('Writing to document');
AOI_Namer(:,2)=[]; % Clears junk data from the second column (actually it kinda just deletes the second columnn)
AOI_Identity(:,2)=[];
%Rewrite it to new doc

%This is all the stuff on sheet 1
xlswrite('Final_Vid', header_1, 'Sheet1', 'A1');
xlswrite('Final_Vid', Start_time_1, 'Sheet1', 'B2');
xlswrite('Final_Vid', hori_posse_1, 'Sheet1', 'C2'); 
xlswrite('Final_Vid', vert_posse_1, 'Sheet1', 'D2');
xlswrite('Final_Vid', fix_frames_1, 'Sheet1', 'E2');
xlswrite('Final_Vid', St1_time, 'Sheet1', 'F2');
xlswrite('Final_Vid', TIME, 'Sheet1', 'G2'); 
xlswrite('Final_Vid', Compare, 'Sheet1', 'H2');
xlswrite('Final_Vid', AOI_Identity, 'Sheet1', 'I2');
xlswrite('Final_Vid', AOI_Namer, 'Sheet1', 'J2');
xlswrite('Final_Vid', generate, 'Sheet1', 'K2');
xlswrite('Final_Vid', Avg_Duartion, 'Sheet1', 'L2');

disp('First Sheet Finished');

%Here begins Sheet 2 (excluding time comparason)
xlswrite('Final_Vid', header_2, 'Sheet2', 'A1');
xlswrite('Final_Vid', Start_time_1, 'Sheet2', 'B2');
xlswrite('Final_Vid', hori_posse_2, 'Sheet2', 'C2');
xlswrite('Final_Vid', vert_posse_2, 'Sheet2', 'D2');
xlswrite('Final_Vid', fix_frames_1, 'Sheet2', 'E2');
xlswrite('Final_Vid', St1_time, 'Sheet2', 'F2');
xlswrite('Final_Vid', TIME, 'Sheet2', 'G2');
xlswrite('Final_Vid', Compare, 'Sheet2', 'H2');
xlswrite('Final_Vid', AOI_Identity, 'Sheet2', 'I2');
xlswrite('Final_Vid', AOI_Namer, 'Sheet2', 'J2');

disp('Second Sheet Completed');

%Here begins Sheet 3 (excluding time comparason)
xlswrite('Final_Vid', header_3, 'Sheet3', 'A1');
xlswrite('Final_Vid', Start_time_1, 'Sheet3', 'B2');
xlswrite('Final_Vid', hori_posse_3, 'Sheet3', 'C2');
xlswrite('Final_Vid', vert_posse_3, 'Sheet3', 'D2');
xlswrite('Final_Vid', fix_frames_1, 'Sheet3', 'E2');
xlswrite('Final_Vid', St1_time, 'Sheet3', 'F2');
xlswrite('Final_Vid', TIME, 'Sheet3', 'G2');
xlswrite('Final_Vid', Compare, 'Sheet3', 'H2');
xlswrite('Final_Vid', AOI_Identity, 'Sheet3', 'I2');
xlswrite('Final_Vid', AOI_Namer, 'Sheet3', 'J2');

disp('Final Sheet Completed');
end