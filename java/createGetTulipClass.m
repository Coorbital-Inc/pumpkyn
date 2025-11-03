function createGetTulipClass()

headerStr = {'public class GetTulip {' 
    '// r = GetTulip.getTulip(6.5, 2, 1);' 
    '//' 
    '// ---- Constants for Earth–Moon CR3BP ----'
    'private static final double muStar = 0.0121505856;'
    'private static final double lStar  = 384400.0;   // km'
    'private static final double tStar  = 2.36059e6;  // s'
    '// ---- Lookup tables ----'
    '// Columns: [tau0, x, z, ydot]'
    ''};


%% Get data for Np = 1:15:
for Np=1:20
    headerStr = genStr(Np,headerStr);
end

%% Write to file:
fid = fopen('getTulipData.txt','w');

for tl=1:size(headerStr,1)
    fprintf(fid,'%s\n',headerStr{tl,1});
end
fclose(fid);

end

function headerStr = genStr(Np,headerStr)
       pm = 1;
[tau0,x0] = pumpkyn.cr3bp.getTulip([],Np,pm);
headerStr{end+1} = ['private static final double[][] tau_x0_N',num2str(Np),' = {'];
for i=1:size(tau0,1)
   if i < size(tau0,1)
        headerStr{end+1,1} = sprintf('{ %.15f, %.15f, %.15f, %.15f },',tau0(i),x0(i,1),x0(i,3),x0(i,5));
   else
        headerStr{end+1,1} = sprintf('{ %.15f, %.15f, %.15f, %.15f }',tau0(i),x0(i,1),x0(i,3),x0(i,5));
   end
end
headerStr{end+1,1} = '};';
headerStr{end+1,1} = '';
end