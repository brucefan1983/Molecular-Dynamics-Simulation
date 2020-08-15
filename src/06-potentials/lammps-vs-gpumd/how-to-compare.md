**比较 GPUMD 和 LAMMPS 中 Tersoff 势函数的力**



​    我们知道，不管用什么 MD 软件，只要用的是同一个势函数和同一组势参数，那么对同样的原子系统，原子的受力情况应该是一样的。否则，至少有一个软件编写错误。最近，我比较了 

GPUMD (https://github.com/brucefan1983/GPUMD) 

中若干势函数和 

LAMMPS (https://github.com/lammps/lammps) 

中对应的势函数，都得到了一致的结果。这里，以硅晶体和 Tersoff 势函数为例做一个记录，防止我自己忘记这个过程。



​    1）写一个 matlab 脚本，产生 GPUMD 和 LAMMPS 分别需要的坐标文件：



clear; close all;

r0 = [0.0, 0.0, 0.5, 0.5, 0.25, 0.25, 0.75, 0.75; ...

   0.0, 0.5, 0.0, 0.5, 0.25, 0.75, 0.25, 0.75; ...

   0.0, 0.5, 0.5, 0.0, 0.25, 0.75, 0.75, 0.25].';

n0 = size(r0, 1);

nxyz = 10 * [1, 1, 1];

N = nxyz(1) * nxyz(2) * nxyz(3) * n0;

a = 5.43 * [1, 1, 1];

box_length = a .* nxyz;



r = zeros(N, 3);

n = 0;

for nx = 0 : nxyz(1) - 1

  for ny = 0 : nxyz(2) - 1

​    for nz = 0 : nxyz(3) - 1

​      for m = 1 : n0

​        n = n + 1;

​        r(n, :) = a .* ([nx,ny,nz] + r0(m, :)) + (rand(1,3)-0.5)*0.1;  % 初始位置有一定的随机性

​      end

​    end

  end

end



% 产生 GPUMD 需要的 xyz.in 文件

fid = fopen('xyz.in', 'w');

fprintf(fid, '%d %d %g\n', N, 4, 3.0);

fprintf(fid, '%d %d %d %g %g %g\n', 1, 1, 1, box_length);

for n =1 : N

  fprintf(fid, '%d %d %g %g %g %g\n', 0, 0, 28, r(n, :));

end

fclose(fid);



% 产生 LAMMPS 需要的 xyz.data 文件

fid = fopen('xyz.data', 'w');

fprintf(fid, '\n%d atoms\n', N);

fprintf(fid, '%d atom types\n', 1);

fprintf(fid, '0 %g xlo xhi\n', box_length(1));

fprintf(fid, '0 %g ylo yhi\n', box_length(2));

fprintf(fid, '0 %g zlo zhi\n', box_length(3));

fprintf(fid, '\nMasses\n\n');

fprintf(fid, '1 28\n');

fprintf(fid, '\nAtoms\n\n');

for n =1 : N

  fprintf(fid, '%d %d %g %g %g\n', n, 1, r(n, :));

end

fclose(fid);



​    2）用 LAMMPS 运行如下脚本，计算一组力（结果保存在 f_lammps.txt 文件中）：



units     metal

read_data xyz.data

pair_style tersoff

\# 文件 SiCGe.tersoff 可在此 (https://github.com/lammps/lammps/tree/master/potentials) 找到

pair_coeff * * SiCGe.tersoff Si(D) 

dump    force all custom 1 f_lammps.txt id fx fy fz

run      0



​    3）用 GPUMD 运行如下脚本，计算另一组力（注意：编译 GPUMD 时，需要加上 -DFORCE 选项）：



potential  potentials/tersoff/si_tersoff_1989.txt  #这个势文件包含在GPUMD的发行版中

\# 用GPUMD计算初始力时只需要一行命令，结果将保存在对应文件夹的 f.out 文件



   4）用如下matlab脚本比较两组力：



clear; close all;

load f.out;

load f_lammps.txt; # 我删掉了开头几行没用的文字（我讨厌这种混合文字与数据的输出文件）

[tmp,index]=sort(f_lammps(:,1));

f_lammps=f_lammps(index,2:4);



figure;

plot(f(:,:)-f_lammps(:,:));

xlabel('Atom Index');

ylabel('F^{gpumd} - F^{lammps}');





我得到的图怎么也传不上来。算了，感兴趣的读者可自己动手尝试。