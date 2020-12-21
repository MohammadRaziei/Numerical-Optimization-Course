function create_report(filename, Title, newton_powel_counter_f, newton_powel_counter_gf, newton_powel_counter_Hf, newton_powel_min_f, newton_powel_min_iter, newton_powel_min_x, newton_rosenbrock_counter_f, newton_rosenbrock_counter_gf, newton_rosenbrock_counter_Hf, newton_rosenbrock_min_f, newton_rosenbrock_min_iter, newton_rosenbrock_min_x, sd_powel_counter_f, sd_powel_counter_gf, sd_powel_min_f, sd_powel_min_iter, sd_powel_min_x, sd_rosenbrock_counter_f, sd_rosenbrock_counter_gf, sd_rosenbrock_min_f, sd_rosenbrock_min_iter, sd_rosenbrock_min_x)
s = @(a)"["+strjoin(string(a(:)),', ')+"]";

readme = [
"# "+Title+"  ";
"  ";
"  ";
"## SD  ";
"  ";
"|   | final x  | final f | # func iter | # func eval | # grad eval | # Hess eval  |";
"|:---:|:---:|:---:|:---:|:---:|:---:|:---:|";
"| Powel | "+s(sd_powel_min_x)+" | "+sd_powel_min_f+" | "+sd_powel_min_iter+" | "+sd_powel_counter_f+" | "+sd_powel_counter_gf+" | - |";
"| Rosenbrock | "+s(sd_rosenbrock_min_x)+" | "+sd_rosenbrock_min_f+" | "+sd_rosenbrock_min_iter+" | "+sd_rosenbrock_counter_f+" | "+sd_rosenbrock_counter_gf+" | - |";
"  ";
"## Newton  ";
"  ";
"|   | final x  | final f | # func iter | # func eval | # grad eval | # Hess eval  |";
"|:---:|:---:|:---:|:---:|:---:|:---:|:---:|";
"| Powel | "+s(newton_powel_min_x)+" | "+newton_powel_min_f+" | "+newton_powel_min_iter+" | "+newton_powel_counter_f+" | "+newton_powel_counter_gf+" | "+newton_powel_counter_Hf+" |";
"| Rosenbrock | "+s(newton_rosenbrock_min_x)+" |"+newton_rosenbrock_min_f+" | "+newton_rosenbrock_min_iter+" | "+newton_rosenbrock_counter_f+" | "+newton_rosenbrock_counter_gf+" | "+newton_rosenbrock_counter_Hf+" |";
];


fileID = fopen(filename,'w');
fprintf(fileID,'%s\n',readme);
fclose(fileID);


end