$ set def $public:[gensimul.example2]
$ run $public:[gensimul]slink
   25432
   20
   2
   0.000
$ copy datafile.msim datafile.dat
$ unknown
$ run $public:[gensimul]msim
$ copy datafile.isim datafile.dat
$ run $public:[gensimul]isim
$ copy datafile.lsim datafile.dat
$ run $public:[gensimul]lsim
