# Binary-Crocodile
Code for poly malware clustering

General purpose of this code is to clustering malware binaries, which were created with
polymorph engines. It’s not guaranteed, that algo will help you with all instances of
specific malware family, but it can help with samples in same generation and with similar
pre-config set. 

The main idea is quite simple - collect entropy characteristic of executable section and map
them to euclid coord system. 

Then map entropy features to harmonic functions characteristics and assemble them with 
fast fourier transform (FFT). It’s not necessary to use FFT. I choose it due to simple 
and flexible mechanic, but you also can choose any math model you like. 

Executable section of malware file contain regions with different entropy - 
native code, packet data, user code, data and so on. 
So, we have three characteristics - size of region, entropy and it offset. 
Size of region == period of harmonic function. 
Offset inside file == phase. 
Value of entropy == amplitude. 

### Entropy tree
![alt tag](https://cloud.githubusercontent.com/assets/1109634/16461894/72c214cc-3e37-11e6-9cba-c541f5a661ee.PNG)

For localization of regions we use binary tree. Initialize root with whole executable sections
then calc entropy of file from start to half for left branch and from half to end. Then repeat
this for both branches while difference of entropy between children's became less than some
fixed value e. I use 0.1, but you can make it even less to reach more strict borders of section. It’s also all about performance. 

![alt tag](https://cloud.githubusercontent.com/assets/1109634/16453137/759655d2-3e14-11e6-9702-5976b0961e81.PNG)

After split, we take crones from each depth. For amplitude of each even region we multiply it
by -1 to provide more effective mapping. For testing I use ~500 malware samples of Backdoor.Bedep for win32. You can find list of
them at file "Win32Bedep". For comparison effectiveness of algo was used simple clustering
based on fuzzy hashing (FH) + Levenshtein distance (LD). I use fixed LD, but is probably
more useful to vary limit based on length of FH. Here is statistic of clustering with FH and entropy maps of groups. 


As you can see, percent of clustered files is very low, but groups is “good”. Files are similar and algo can “ignore” minor
differences in file structure. Here is example of mapping entropy features for few clusterized with FH+LD.

### Percentage of clustered with FH+LD files.
![alt tag](https://cloud.githubusercontent.com/assets/1109634/16453144/7c21a334-3e14-11e6-8b36-7afa0bd1d08a.png)


### Clusters FH+LD
![alt tag](https://cloud.githubusercontent.com/assets/1109634/16453106/5ceb6a7c-3e14-11e6-95a8-bd811408c448.png)
![alt tag](https://cloud.githubusercontent.com/assets/1109634/16453108/6083acf8-3e14-11e6-9f1a-2ad6828aba52.png)

### Raw entropy data, without mapping.
![alt tag](https://cloud.githubusercontent.com/assets/1109634/16461756/d4c869ba-3e36-11e6-918b-48da9371b6f8.png)
![alt tag](https://cloud.githubusercontent.com/assets/1109634/16461757/d4e1f9c0-3e36-11e6-975e-6c6d61256646.png)


For comparison of “entropy” function we integrate them and compare scales with fixed limit (0.05).

### Fourie functions for groups

![alt tag](https://cloud.githubusercontent.com/assets/1109634/16462183/c485543a-3e38-11e6-9c5e-9476eb2c66ac.png)
![alt tag](https://cloud.githubusercontent.com/assets/1109634/16462182/c483f00e-3e38-11e6-87d1-4d841f0b9692.png)
![alt tag](https://cloud.githubusercontent.com/assets/1109634/16462181/c46e15ae-3e38-11e6-9b56-9274117fd89d.png)

### Entropy maps for groups

![alt tag](https://cloud.githubusercontent.com/assets/1109634/16462130/95ea3e42-3e38-11e6-9943-5b5b8d033cd1.png)
![alt tag](https://cloud.githubusercontent.com/assets/1109634/16462131/95f053ae-3e38-11e6-915b-e553bd04f5ca.png)
