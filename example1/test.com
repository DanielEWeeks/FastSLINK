$ set def $public:[gensimul.example1]
$ run $public:[gensimul]slink
   29876
   200
   1
   0.000
$ copy datafile.msim datafile.dat
$ unknown
$ run $public:[gensimul]msim
$ copy datafile.isim datafile.dat
$ run $public:[gensimul]isim
$ copy datafile.lsim datafile.dat
$ run $public:[gensimul]lsim
