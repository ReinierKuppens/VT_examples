
function mat = str2inc(str)

mS  = max(str);
M   = incPattern2(mS);
mat = M(:,str);