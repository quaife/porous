for j = 1:470
  fileName = ['inner' num2str(j,'%03d') '.dat'];

  fid = fopen(fileName,'w');
  for k=1:128 
    str=[num2str(x(k,j),'%8.4e') ' & ' num2str(y(k,j),'%8.4e')]; 
    fprintf(fid,'%s\n',str); 
  end;
  str=[num2str(x(1,j),'%8.4e') ' & ' num2str(y(1,j),'%8.4e')];
  fprintf(fid,'%s\n',str);

end

