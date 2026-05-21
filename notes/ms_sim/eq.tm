<TeXmacs|2.1.4>

<style|<tuple|generic|chinese>>

<\body>
  <\hide-preamble>
    \;

    <assign||<\macro>
      <with|font-family|rm| d>
    </macro>>

    <assign|ud|<\macro>
      <with|font-family|rm| d>

      \;
    </macro>>

    <assign|edit-math|<macro|body|<math|<arg|body>>>>

    <assign|ui|<\macro>
      <with|font-family|rm|i>
    </macro>>

    <assign|ue|<\macro>
      <with|font-family|rm|e>
    </macro>>
  </hide-preamble>

  Area term when the two ions experience different amplitude modulations

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<gamma\><rsup|n><rsub|k>>|<cell|\<equiv\>
    >|<cell| <big|int><rsup|t<rsub|n+1>><rsub|t<rsub|n>>dt
    \<Omega\><rsub|1><around*|(|t|)>e<rsup|i
    \<theta\><rsub|k><around*|(|t|)>><big|int><rsub|t<rsub|n>><rsup|t>dt<rprime|'>
    \<Omega\><rsub|2><around*|(|t<rprime|'>|)>e<rsup|-i\<theta\><rsub|k><around*|(|t<rprime|'>|)>>>>|<row|<cell|>|<cell|=>|<cell|
    <big|int><rsup|\<tau\><rsub|n>><rsub|0>dt
    <around*|(|\<Omega\><rsub|1n>+\<Omega\><rsub|1n><rprime|'>t|)>e<rsup|i
    \<delta\><rsup|n><rsub|k>t><big|int><rsub|0><rsup|t>dt<rprime|'>
    <around*|(|\<Omega\><rsub|2n>+\<Omega\><rsub|2n><rprime|'>t<rprime|'>|)>e<rsup|-i
    \<delta\><rsup|n><rsub|k>t<rprime|'>>>>|<row|<cell|>|<cell|=>|<cell|
    <big|int><rsup|\<tau\><rsub|n>><rsub|0>dt
    <around*|(|\<Omega\><rsub|1n>+\<Omega\><rsub|1n><rprime|'>t|)>e<rsup|i
    \<delta\><rsup|n><rsub|k>t><around*|(|<around*|(|<frac|
    \<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsub|k><rsup|n>|\<delta\><rsub|k><rsup|n>
    <rsup|2>>+i<frac| \<Omega\><rsub|2n><rprime|'>t|\<delta\><rsub|k><rsup|n>>|)>e<rsup|-i\<delta\><rsup|n><rsub|k>t>-<frac|
    \<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsup|n><rsub|k>|\<delta\><rsup|n><rsub|k>
    <rsup|2>>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<delta\><rsup|n><rsub|k>
    <rsup|2>> <big|int><rsup|\<tau\><rsub|n>><rsub|0>dt
    \<Omega\><rsub|1n><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsub|k><rsup|n>|)>+<around*|(|\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n>\<delta\><rsub|k><rsup|n>+i\<Omega\><rsub|1n>\<Omega\><rsub|2n><rprime|'>\<delta\><rsub|k><rsup|n>|)>t+\<Omega\><rsub|1n><rprime|'>
    i\<Omega\><rsub|2n><rprime|'>\<delta\><rsub|k><rsup|n>t<rsup|2>-<around*|(|\<Omega\><rsub|1n>+\<Omega\><rsub|1n><rprime|'>t|)><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsup|n><rsub|k>|)>e<rsup|i
    \<delta\><rsup|n><rsub|k>t>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<delta\><rsup|n><rsub|k>
    <rsup|2>> <around*|(|\<Omega\><rsub|1n><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsub|k><rsup|n>|)>\<tau\><rsub|n>+<frac|1|2><around*|(|\<Omega\><rsub|1n><rprime|'>\<Omega\>
    <rsub|2n><rprime|'>+i\<Omega\><rsub|1n><rprime|'>\<Omega\>
    <rsub|2n>\<delta\><rsub|k><rsup|n>+i\<Omega\><rsub|1n>\<Omega\><rsub|2n><rprime|'>\<delta\><rsub|k><rsup|n>|)>\<tau\><rsup|2><rsub|n<rsup|>>+<frac|i|3>\<Omega\><rsub|1n><rprime|'>
    \<Omega\><rsub|2n><rprime|'>\<delta\><rsub|k><rsup|n>\<tau\><rsub|n><rsup|3>-\<Omega\><rsub|1n><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsup|n><rsub|k>|)><around*|(|F<rsub|0><around*|(|\<delta\><rsub|k><rsup|n>,\<tau\><rsub|n>|)>-F<rsub|0><around*|(|\<delta\><rsub|k><rsup|n>,0|)>|)>-\<Omega\><rsub|1n><rprime|'><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsup|n><rsub|k>|)><around*|(|F<rsub|1><around*|(|\<delta\><rsub|k><rsup|n>,\<tau\><rsub|n>|)>-F<rsub|1><around*|(|\<delta\><rsub|k><rsup|n>,0|)>|)>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<delta\><rsup|n><rsub|k>
    <rsup|2>> <around*|(|\<Omega\><rsub|1n><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsub|k><rsup|n>|)>\<tau\><rsub|n>+<frac|1|2><around*|(|\<Omega\><rsub|1n><rprime|'>\<Omega\>
    <rsub|2n><rprime|'>+i\<Omega\><rsub|1n><rprime|'>\<Omega\>
    <rsub|2n>\<delta\><rsub|k><rsup|n>+i\<Omega\><rsub|1n>\<Omega\><rsub|2n><rprime|'>\<delta\><rsub|k><rsup|n>|)>\<tau\><rsup|2><rsub|n<rsup|>>+<frac|i|3>\<Omega\><rsub|1n><rprime|'>
    \<Omega\><rsub|2n><rprime|'> \<delta\><rsub|k><rsup|n>\<tau\><rsub|n><rsup|3>+<frac|i|\<delta\><rsub|k><rsup|n>>
    \<Omega\><rsub|1n><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsup|n><rsub|k>|)><around*|(|e<rsup|i\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>>-1|)>-\<Omega\><rsub|1n><rprime|'><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsup|n><rsub|k>|)><around*|(|<frac|-i\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>+1|\<delta\><rsub|k><rsup|n>
    <rsup|2>>e<rsup|i\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>>-<frac|1|\<delta\><rsub|k><rsup|n>
    <rsup|2>>|)>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|\<delta\><rsup|n><rsub|k>
    <rsup|2>> <around*|(|\<Omega\><rsub|1n><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsub|k><rsup|n>|)>\<tau\><rsub|n>+<frac|1|2><around*|(|\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n>\<delta\><rsub|k><rsup|n>+i\<Omega\><rsub|1n>\<Omega\><rsub|2n><rprime|'>\<delta\><rsub|k><rsup|n>|)>\<tau\><rsup|2><rsub|n<rsup|>>+<frac|i|3>\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n><rprime|'>\<delta\><rsub|k><rsup|n>\<tau\><rsub|n><rsup|3>-<frac|1|\<delta\><rsub|k><rsup|n>
    <rsup|2>> <around*|(|\<Omega\><rsub|1n>\<Omega\><rsub|2n>\<delta\><rsup|n><rsub|k>
    <rsup|2>+\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n><rprime|'>|)><around*|(|e<rsup|i\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>>-1|)>+<frac|i|\<delta\><rsub|k><rsup|n>>
    <around*|(|\<Omega\><rsub|1n>\<Omega\><rsub|2n><rprime|'>-\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n>|)><around*|(|e<rsup|i\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>>-1|)>+i
    \<Omega\><rsub|1n><rprime|'><around*|(|\<Omega\><rsub|2n><rprime|'>+i\<Omega\><rsub|2n>\<delta\><rsup|n><rsub|k>|)><frac|\<tau\><rsub|n>|\<delta\><rsub|k><rsup|n>>e<rsup|i\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>>|)>>>|<row|<cell|>|<cell|=>|<cell|\<Omega\><rsub|1n>\<Omega\><rsub|2n>\<tau\><rsub|n><rsup|2><frac|1-cos<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>|<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>
    <rsup|2>>+\<Omega\><rsub|1n>\<Omega\><rsub|2n><rprime|'>\<tau\><rsub|n><rsup|3><frac|\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>-sin<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>|)><rsup|3>>>>|<row|<cell|>|<cell|>|<cell|+\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n>\<tau\><rsub|n><rsup|3><frac|
    sin<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>-\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>cos<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>|)><rsup|3>>>>|<row|<cell|>|<cell|>|<cell|+\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n><rprime|'>\<tau\><rsup|4><rsub|n<rsup|>><frac|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsup|><rsub|n>|)><rsup|2>/2+1-cos<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>-\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>sin<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsup|><rsub|n>|)><rsup|4>>>>|<row|<cell|>|<cell|>|<cell|+i\<Omega\><rsub|1n>\<Omega\><rsub|2n>\<tau\><rsub|n><rsup|2><frac|\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>-sin<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>|<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>
    <rsup|2>>+i\<Omega\><rsub|1n>\<Omega\><rsub|2n><rprime|'>\<tau\><rsub|n><rsup|3><frac|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsup|><rsub|n<rsup|>>|)><rsup|2>/2+cos<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>-1|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsup|><rsub|n<rsup|>>|)><rsup|3>>>>|<row|<cell|>|<cell|>|<cell|+i
    \<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n>\<tau\><rsub|n><rsup|3><frac|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsup|><rsub|n<rsup|>>|)><rsup|2>/2+1-cos<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>-\<delta\><rsub|k><rsup|n>\<tau\><rsup|><rsub|n<rsup|>>sin<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsup|><rsub|n<rsup|>>|)><rsup|3>>
    >>|<row|<cell|>|<cell|>|<cell|+i\<Omega\><rsub|1n><rprime|'>\<Omega\><rsub|2n><rprime|'>\<tau\><rsub|n><rsup|4>
    <frac|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>|)><rsup|3>/3-sin<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>+\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>cos<around*|(|\<delta\><rsup|n><rsub|k>\<tau\><rsub|n>|)>|<around*|(|\<delta\><rsub|k><rsup|n>\<tau\><rsub|n>|)><rsup|4>>>>>>
  </eqnarray*>

  \;

  With different frequencies as well,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<gamma\><rsup|n><rsub|k>>|<cell|\<equiv\>
    >|<cell| <big|int><rsup|t<rsub|n+1>><rsub|t<rsub|n>>dt \<Omega\>
    <rsub|1><around*|(|t|)>e<rsup|i \<theta\><rsub|1k><around*|(|t|)>><big|int><rsub|t<rsub|n>><rsup|t>dt<rprime|'>
    \<Omega\><rsub|2><around*|(|t<rprime|'>|)>e<rsup|-i\<theta\><rsub|2k><around*|(|t<rprime|'>|)>>>>|<row|<cell|>|<cell|=>|<cell|
    <big|int><rsup|\<tau\><rsub|n>><rsub|0>dt
    <around*|(|\<Omega\><rsub|1n>+\<Omega\><rsub|1n><rprime|'>t|)>e<rsup|i
    \<delta\><rsup|1n><rsub|k>t><big|int><rsub|0><rsup|t>dt<rprime|'>
    <around*|(|\<Omega\><rsub|2n>+\<Omega\><rsub|2n><rprime|'>t<rprime|'>|)>e<rsup|-i
    \<delta\><rsup|2n><rsub|k>t<rprime|'>>>>>>
  </eqnarray*>

  Simplifying the symbols,

  <\eqnarray>
    <tformat|<table|<row|<cell|d<rsub|<around*|{|1,2|}>>>|<cell|\<equiv\>>|<cell|\<delta\><rsup|<around*|{|1,2|}>n><rsub|k>\<tau\><rsub|n>>>|<row|<cell|o<rsub|<around*|{|1,2|}>>>|<cell|\<equiv\>>|<cell|\<Omega\><rsub|<around*|{|1,2|}>n>\<tau\><rsub|n>>>|<row|<cell|o<rprime|'><rsub|<around*|{|1,2|}>>>|<cell|\<equiv\>>|<cell|\<Omega\><rsub|<around*|{|1,2|}>n><rprime|'>\<tau\><rsub|n><rsup|2>>>>>
  </eqnarray>

  <\eqnarray>
    <tformat|<table|<row|<cell|\<gamma\><rsup|n><rsub|k>>|<cell|=>|<cell|
    <big|int><rsup|1><rsub|0>dt \ <around*|(|\<Omega\><rsub|1n>\<tau\><rsub|n>+\<Omega\><rsub|1n><rprime|'>\<tau\><rsup|2><rsub|n>t|)>e<rsup|i
    \<delta\><rsup|1n><rsub|k>\<tau\><rsub|n>t><big|int><rsub|0><rsup|t>dt<rprime|'>
    <around*|(|\<Omega\><rsub|2n>\<tau\><rsub|n>+\<Omega\><rsub|2n><rprime|'>\<tau\><rsub|n><rsup|2>t<rprime|'>|)>e<rsup|-i
    \<delta\><rsup|2n><rsub|k>\<tau\><rsub|n>t<rprime|'>>>>|<row|<cell|>|<cell|=>|<cell|<big|int><rsup|1><rsub|0>dt
    \ <around*|(|o<rsub|1>+o<rsub|1><rprime|'>t|)>e<rsup|i
    d<rsub|1>t><big|int><rsub|0><rsup|t>dt<rprime|'>
    <around*|(|o<rsub|2>+o<rprime|'><rsub|2>t<rprime|'>|)>e<rsup|-i
    d<rsub|2>t<rprime|'>>>>|<row|<cell|>|<cell|=>|<cell|<big|int><rsup|1><rsub|0>dt
    \ <around*|(|o<rsub|1>+o<rsub|1><rprime|'>t|)>e<rsup|i
    d<rsub|1>t><around*|(|<frac|i o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>+i
    o<rprime|'><rsub|2>d<rsub|2> t|d<rsub|2<rsup|>><rsup|2>>e<rsup|-i
    d<rsub|2>t>-<frac|i o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|d<rsub|2<rsup|>><rsup|2>>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|d<rsub|2><rsup|2>><big|int><rsup|1><rsub|0>dt
    \ <around*|(|o<rsub|1>+o<rsub|1><rprime|'>t|)><around*|(|i
    o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>+i o<rprime|'><rsub|2>d<rsub|2>
    t|)>e<rsup|i<around*|(| d<rsub|2>-d<rsub|1>|)>t>-<around*|(|o<rsub|1>+o<rsub|1><rprime|'>t|)><around*|(|i
    o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)>e<rsup|i
    d<rsub|1>t>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|d<rsub|2><rsup|2>><big|int><rsup|1><rsub|0>dt
    o<rsub|1><around*|(|i o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)>e<rsup|i<around*|(|
    d<rsub|2>-d<rsub|1>|)>t>+<around*|(|o<rsub|1><rprime|'>o<rprime|'><rsub|2>+i
    o<rsub|1><rprime|'>o<rsub|2>d<rsub|2>+i
    o<rsub|1>o<rprime|'><rsub|2>d<rsub|2>|)>t e<rsup|i<around*|(|
    d<rsub|2>-d<rsub|1>|)>t>+i o<rsub|1><rprime|'>o<rprime|'><rsub|2>d<rsub|2>t<rsup|2>e<rsup|i<around*|(|
    d<rsub|2>-d<rsub|1>|)>t>-o<rsub|1><around*|(|i
    o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)>e<rsup|i
    d<rsub|1>t>-o<rsub|1><rprime|'><around*|(|i
    o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)>t e<rsup|i
    d<rsub|1>t>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|d<rsub|2><rsup|2>><around*|(|i<rsub|>o<rsub|1><around*|(|i<rsub|>o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)><frac|1-e<rsup|i<around*|(|d<rsub|2>-d<rsub|1>|)>>|d<rsub|2>-d<rsub|1>>+<around*|(|o<rsub|1><rprime|'>o<rprime|'><rsub|2>+i
    o<rsub|1><rprime|'>o<rsub|2>d<rsub|2>+i
    o<rsub|1>o<rprime|'><rsub|2>d<rsub|2>|)><frac|-1+<around*|(|1-i<around*|(|d<rsub|2>-d<rsub|1>|)>|)>e<rsup|i<around*|(|d<rsub|2>-d<rsub|1>|)>>|<around*|(|d<rsub|2>-d<rsub|1>|)><rsup|2>>+i
    o<rsub|1><rprime|'>o<rprime|'><rsub|2>d<rsub|2><frac|<around*|(|2i+2<around*|(|d<rsub|2>-d<rsub|1>|)>-i<around*|(|d<rsub|2>-d<rsub|1>|)><rsup|2>|)>e<rsup|i<around*|(|d<rsub|2>-d<rsub|1>|)>>-2i|<around*|(|d<rsub|2>-d<rsub|1>|)><rsup|3>>-i
    o<rsub|1><around*|(|i o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)><frac|1-e<rsup|i
    d<rsub|1>>|d<rsub|1>>-o<rsub|1><rprime|'><around*|(|i
    o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)><frac|-1+<around*|(|1-i
    d<rsub|1>|)>e<rsup|i d<rsub|1>>|d<rsub|1><rsup|2>>|)>>>>>
  </eqnarray>

  Define <math|\<Delta\>\<equiv\>d<rsub|2>-d<rsub|1>>

  <\eqnarray>
    <tformat|<table|<row|<cell|\<gamma\><rsup|n><rsub|k>>|<cell|=>|<cell|<frac|1|d<rsub|2><rsup|2>><around*|(|i<rsub|>o<rsub|1><around*|(|i<rsub|>o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)><frac|1-e<rsup|i\<Delta\>>|\<Delta\>>+<around*|(|o<rsub|1><rprime|'>o<rprime|'><rsub|2>+i<rsub|>o<rsub|1><rprime|'>o<rsub|2>d<rsub|2>+i<rsub|>o<rsub|1>o<rprime|'><rsub|2>d<rsub|2>|)><frac|-1+<around*|(|1-i\<Delta\>|)>e<rsup|i\<Delta\>>|\<Delta\><rsup|2>>+i
    o<rsub|1><rprime|'>o<rprime|'><rsub|2>d<rsub|2><frac|<around*|(|2i+2\<Delta\>-i\<Delta\><rsup|2>|)>e<rsup|i\<Delta\>>-2i|\<Delta\><rsup|3>>-i
    o<rsub|1><around*|(|i<rsub|>o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)><frac|1-e<rsup|i
    d<rsub|1>>|d<rsub|1>>-o<rsub|1><rprime|'><around*|(|i
    o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)><frac|-1+<around*|(|1-i
    d<rsub|1>|)>e<rsup|i d<rsub|1>>|d<rsub|1><rsup|2>>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|1|d<rsub|2><rsup|2>><around*|(|i<rsub|>o<rsub|1><around*|(|i<rsub|>o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)><frac|1-cos\<Delta\>-i
    sin\<Delta\>|\<Delta\>>+<around*|(|o<rsub|1><rprime|'>o<rprime|'><rsub|2>+i
    o<rsub|1><rprime|'>o<rsub|2>d<rsub|2>+i
    o<rsub|1>o<rprime|'><rsub|2>d<rsub|2>|)><frac|-1+<around*|(|1-i\<Delta\>|)><around*|(|cos\<Delta\>+i
    sin\<Delta\>|)>|\<Delta\><rsup|2>>+o<rsub|1><rprime|'>o<rprime|'><rsub|2>d<rsub|2><frac|<around*|(|-2+2i\<Delta\>+\<Delta\><rsup|2>|)><around*|(|cos\<Delta\>+i
    sin\<Delta\>|)>+2|\<Delta\><rsup|3>>+o<rsub|1><around*|(|o<rsub|2>d<rsub|2>-i
    o<rprime|'><rsub|2>|)><frac|1-cos<around*|(|d<rsub|1 >|)>-i
    sin<around*|(|d<rsub|1>|)>|d<rsub|1>>-o<rsub|1><rprime|'><around*|(|i
    o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)><frac|-1+<around*|(|1-i
    d<rsub|1>|)><around*|(|cos<around*|(|d<rsub|1>|)>+i
    sin<around*|(|d<rsub|1>|)>|)>|d<rsub|1><rsup|2>>|)>>>|<row|<cell|>|<cell|=>|<cell|-o<rsub|1>o<rsub|2><frac|1-cos\<Delta\>-i
    sin\<Delta\>|\<Delta\>d<rsub|2>>+i<rsub|>o<rsub|1>o<rprime|'><rsub|2><frac|1-cos\<Delta\>-i
    sin\<Delta\>|\<Delta\>d<rsub|2><rsup|2>>+<around*|(|o<rsub|1><rprime|'>o<rprime|'><rsub|2>+i
    o<rsub|1><rprime|'>o<rsub|2>d<rsub|2>+i
    o<rsub|1>o<rprime|'><rsub|2>d<rsub|2>|)><frac|cos\<Delta\>+\<Delta\>sin\<Delta\>-1+i
    sin\<Delta\>-i\<Delta\>cos\<Delta\>|\<Delta\><rsup|2>d<rsub|2><rsup|2>>+o<rsub|1><rprime|'>o<rprime|'><rsub|2><frac|<around*|(|\<Delta\><rsup|2>-2|)>cos\<Delta\>+i<around*|(|\<Delta\><rsup|2>-2|)>sin\<Delta\>+2i\<Delta\>
    cos\<Delta\>- 2\<Delta\>sin\<Delta\>+2|\<Delta\><rsup|3>d<rsub|2>>+o<rsub|1><around*|(|o<rsub|2>d<rsub|2>-i
    o<rprime|'><rsub|2>|)><frac|1-cos<around*|(|d<rsub|1 >|)>-i
    sin<around*|(|d<rsub|1>|)>|d<rsub|1>d<rsub|2><rsup|2>>-o<rsub|1><rprime|'><around*|(|i
    o<rsub|2>d<rsub|2>+o<rprime|'><rsub|2>|)><frac|-1+<around*|(|1-i
    d<rsub|1>|)><around*|(|cos<around*|(|d<rsub|1>|)>+i
    sin<around*|(|d<rsub|1>|)>|)>|d<rsub|1><rsup|2>d<rsub|2><rsup|2>>>>|<row|<cell|>|<cell|=>|<cell|<frac|o<rsub|1>o<rsub|2>|d<rsub|2>><around*|(|<frac|1-cos<around*|(|d<rsub|1
    >|)>-i sin<around*|(|d<rsub|1>|)>|d<rsub|1>>-<frac|1-cos\<Delta\>-i
    sin\<Delta\>|\<Delta\>>|)>>>|<row|<cell|>|<cell|>|<cell|+<frac|o<rsub|1>o<rprime|'><rsub|2>|d<rsub|2>><around*|(|<frac|i
    cos\<Delta\>+i\<Delta\>sin\<Delta\>-i-
    sin\<Delta\>+\<Delta\>cos\<Delta\>|\<Delta\><rsup|2>>-<frac|i-i
    cos<around*|(|d<rsub|1 >|)>+ sin<around*|(|d<rsub|1>|)>|d<rsub|1>d<rsub|2>>+<rsub|><frac|i-i<rsub|>cos\<Delta\>+sin\<Delta\>|\<Delta\>d<rsub|2>>|)>>>|<row|<cell|>|<cell|>|<cell|+<frac|o<rsub|1><rprime|'>o<rsub|2>|d<rsub|2>><around*|(|<frac|i
    cos\<Delta\>+i\<Delta\>sin\<Delta\>-i-
    sin\<Delta\>+\<Delta\>cos\<Delta\>|\<Delta\><rsup|2>>-<frac|-i+i
    cos<around*|(|d<rsub|1>|)>+ d<rsub|1>cos<around*|(|d<rsub|1>|)>-sin<around*|(|d<rsub|1>|)>+
    i d<rsub|1>sin<around*|(|d<rsub|1>|)>|d<rsub|1><rsup|2>>|)>>>|<row|<cell|>|<cell|>|<cell|+<frac|o<rsub|1><rprime|'>o<rprime|'><rsub|2>|d<rsub|2>><around*|(|<frac|cos\<Delta\>+\<Delta\>sin\<Delta\>-1+i
    sin\<Delta\>-i\<Delta\>cos\<Delta\>|\<Delta\><rsup|2>d<rsub|2>>+<frac|<around*|(|\<Delta\><rsup|2>-2|)>cos\<Delta\>+i<around*|(|\<Delta\><rsup|2>-2|)>sin\<Delta\>+2i\<Delta\>
    cos\<Delta\>- 2\<Delta\>sin\<Delta\>+2|\<Delta\><rsup|3>>-<frac|-1+cos<around*|(|d<rsub|1>|)>-i
    d<rsub|1>cos<around*|(|d<rsub|1>|)>+i sin<around*|(|d<rsub|1>|)>+
    d<rsub|1>sin<around*|(|d<rsub|1>|)>|d<rsub|1><rsup|2>d<rsub|2>>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|o<rsub|1>o<rsub|2>|d<rsub|2>><around*|(|<frac|1-cos<around*|(|d<rsub|1
    >|)>|d<rsub|1>>-<frac|1-cos\<Delta\>|\<Delta\>>|)>>>|<row|<cell|>|<cell|>|<cell|+<frac|o<rsub|1>o<rprime|'><rsub|2>|d<rsub|2>><around*|(|<frac|\<Delta\>cos\<Delta\>-
    sin\<Delta\>|\<Delta\><rsup|2>>-<frac|sin<around*|(|d<rsub|1>|)>|d<rsub|1>d<rsub|2>>+<rsub|><frac|sin\<Delta\>|\<Delta\>d<rsub|2>>|)>>>|<row|<cell|>|<cell|>|<cell|+<frac|o<rsub|1><rprime|'>o<rsub|2>|d<rsub|2>><around*|(|<frac|\<Delta\>cos\<Delta\>-
    sin\<Delta\>|\<Delta\><rsup|2>>-<frac|d<rsub|1>cos<around*|(|d<rsub|1>|)>-sin<around*|(|d<rsub|1>|)>|d<rsub|1><rsup|2>>|)>>>|<row|<cell|>|<cell|>|<cell|+<frac|o<rsub|1><rprime|'>o<rprime|'><rsub|2>|d<rsub|2>><around*|(|<frac|cos\<Delta\>+\<Delta\>sin\<Delta\>-1|\<Delta\><rsup|2>d<rsub|2>>+<frac|<around*|(|\<Delta\><rsup|2>-2|)>cos\<Delta\>-
    2\<Delta\>sin\<Delta\>+2|\<Delta\><rsup|3>>-<frac|-1+cos<around*|(|d<rsub|1>|)>+
    d<rsub|1>sin<around*|(|d<rsub|1>|)>|d<rsub|1><rsup|2>d<rsub|2>>|)>>>|<row|<cell|>|<cell|>|<cell|+i<frac|o<rsub|1>o<rsub|2>|d<rsub|2>><around*|(|<frac|-sin<around*|(|d<rsub|1>|)>|d<rsub|1>>+<frac|sin\<Delta\>|\<Delta\>>|)>>>|<row|<cell|>|<cell|>|<cell|+i<frac|o<rsub|1>o<rprime|'><rsub|2>|d<rsub|2>><around*|(|<frac|cos\<Delta\>+\<Delta\>sin\<Delta\>-1|\<Delta\><rsup|2>>-<frac|1-cos<around*|(|d<rsub|1
    >|)>|d<rsub|1>d<rsub|2>>+<rsub|><frac|1-<rsub|>cos\<Delta\>|\<Delta\>d<rsub|2>>|)>>>|<row|<cell|>|<cell|>|<cell|+i<frac|o<rsub|1><rprime|'>o<rsub|2>|d<rsub|2>><around*|(|<frac|cos\<Delta\>+\<Delta\>sin\<Delta\>-1|\<Delta\><rsup|2>>-<frac|cos<around*|(|d<rsub|1>|)>+
    d<rsub|1>sin<around*|(|d<rsub|1>|)>-1|d<rsub|1><rsup|2>>|)>>>|<row|<cell|>|<cell|>|<cell|+i<frac|o<rsub|1><rprime|'>o<rprime|'><rsub|2>|d<rsub|2>><around*|(|<frac|sin\<Delta\>-\<Delta\>cos\<Delta\>|\<Delta\><rsup|2>d<rsub|2>>+<frac|<around*|(|\<Delta\><rsup|2>-2|)>sin\<Delta\>+2\<Delta\>
    cos\<Delta\>|\<Delta\><rsup|3>>-<frac|-
    d<rsub|1>cos<around*|(|d<rsub|1>|)>+sin<around*|(|d<rsub|1>|)>|d<rsub|1><rsup|2>d<rsub|2>>|)>>>>>
  </eqnarray>

  Check if <with|mode|math|d<rsub|1>=d<rsub|2>=d>

  <\eqnarray>
    <tformat|<table|<row|<cell|\<gamma\><rsub|k><rsup|n>>|<cell|=>|<cell|o<rsub|1>o<rsub|2><frac|1-cos<around*|(|d<rsub|
    >|)>|d<rsup|2><rsub|>>+o<rsub|1>*o<rprime|'><rsub|2><frac|d-sin<around*|(|d|)>|d<rsup|3><rsub|>>>>|<row|<cell|>|<cell|>|<cell|+o<rsub|1><rprime|'>o<rsub|2><frac|sin<around*|(|d<rsub|>|)>-d<rsub|>cos<around*|(|d<rsub|>|)>|d<rsup|3><rsub|>>>>|<row|<cell|>|<cell|>|<cell|+o<rsub|1><rprime|'>o<rprime|'><rsub|2><frac|d<rsup|2>/2+1-cos<around*|(|d|)>-d*sin<around*|(|d|)>|d<rsup|4>>>>|<row|<cell|>|<cell|>|<cell|+i*o<rsub|1>o<rsub|2><frac|d-sin<around*|(|d|)>|d<rsup|2>>+i<rsub|>o<rsub|1>o<rprime|'><rsub|2><frac|d<rsup|2>/2+cos<around*|(|d|)>-1|d<rsup|3>>>>|<row|<cell|>|<cell|>|<cell|+i*o<rsub|1><rprime|'>o<rsub|2>
    <frac|d<rsup|2>/2+1-cos<around*|(|d<rsub|>|)>-
    d<rsub|>sin<around*|(|d|)>|d<rsup|3>>>>|<row|<cell|>|<cell|>|<cell|+i*o<rsub|1><rprime|'>o<rprime|'><rsub|2><frac|d<rsup|3>/3+
    d<rsub|>cos<around*|(|d|)>-sin<around*|(|d|)>|d<rsup|4>>>>>>
  </eqnarray>

  which agrees with the previous result.
</body>

<\initial>
  <\collection>
    <associate|page-height|auto>
    <associate|page-medium|paper>
    <associate|page-type|letter>
    <associate|page-width|auto>
  </collection>
</initial>