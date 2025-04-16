<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Method of Moments for COGARCH(P,Q) — gmm • The YUIMA Project</title><!-- jquery --><script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js" integrity="sha512-v2CJ7UaYy4JwqLDIrZUI/4hqeoQieOmAZNXBeQyjo21dadnwR+8ZaIJVT8EE2iyI61OV8e6M8PP2/4hpQINQ/g==" crossorigin="anonymous" referrerpolicy="no-referrer"></script><!-- Bootstrap --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha256-bZLfwXAP04zRMK2BjiO8iu9pf4FbLqX6zitd+tIvLhE=" crossorigin="anonymous"><script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha256-nuL8/2cJ5NDSSwnKD8VqreErSWHtnEP9E7AySL+1ev4=" crossorigin="anonymous"></script><!-- bootstrap-toc --><link rel="stylesheet" href="../bootstrap-toc.css"><script src="../bootstrap-toc.js"></script><!-- Font Awesome icons --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/all.min.css" integrity="sha256-mmgLkCYLUQbXn0B1SRqzHar6dCnv9oZFPEC1g1cwlkk=" crossorigin="anonymous"><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/v4-shims.min.css" integrity="sha256-wZjR52fzng1pJHwx4aV2AO3yyTOXrcDW7jBpJtTwVxw=" crossorigin="anonymous"><!-- clipboard.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/clipboard.js/2.0.6/clipboard.min.js" integrity="sha256-inc5kl9MA1hkeYUt+EC3BhlIgyp/2jDIyBLS6k3UxPI=" crossorigin="anonymous"></script><!-- headroom.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/headroom.min.js" integrity="sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/jQuery.headroom.min.js" integrity="sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=" crossorigin="anonymous"></script><!-- pkgdown --><link href="../pkgdown.css" rel="stylesheet"><script src="../pkgdown.js"></script><meta property="og:title" content="Method of Moments for COGARCH(P,Q) — gmm"><meta property="og:description" content="The function returns the estimated parameters of a COGARCH(P,Q) model. The parameters are obtained by matching theoretical vs empirical autocorrelation function. The theoretical autocorrelation function is computed according the methodology developed in Chadraa (2009)."><!-- mathjax --><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" integrity="sha256-nvJJv9wWKEm88qvoQl9ekL2J+k/RWIsaSScxxlsrv8k=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/config/TeX-AMS-MML_HTMLorMML.js" integrity="sha256-84DKXVJXs0/F8OTMzX4UR909+jtl4G7SPypPavF+GfA=" crossorigin="anonymous"></script><!--[if lt IE 9]>
<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
<![endif]--></head><body data-spy="scroll" data-target="#toc">


    <div class="container template-reference-topic">
      <header><div class="navbar navbar-default navbar-fixed-top" role="navigation">
  <div class="container">
    <div class="navbar-header">
      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false">
        <span class="sr-only">Toggle navigation</span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
      </button>
      <span class="navbar-brand">
        <a class="navbar-link" href="../index.html">The YUIMA Project</a>
        <span class="version label label-default" data-toggle="tooltip" data-placement="bottom" title="">1.15.30</span>
      </span>
    </div>

    <div id="navbar" class="navbar-collapse collapse">
      <ul class="nav navbar-nav"><li>
  <a href="../reference/index.html">Reference</a>
</li>
      </ul><ul class="nav navbar-nav navbar-right"><li>
  <a href="https://github.com/yuimaproject/yuima/" class="external-link">
    <span class="fab fa-github fa-lg"></span>

  </a>
</li>
      </ul></div><!--/.nav-collapse -->
  </div><!--/.container -->
</div><!--/.navbar -->



      </header><div class="row">
  <div class="col-md-9 contents">
    <div class="page-header">
    <h1>Method of Moments for COGARCH(P,Q)</h1>

    <div class="hidden name"><code>gmm.rd</code></div>
    </div>

    <div class="ref-description">
    <p>The function returns the estimated parameters of a COGARCH(P,Q) model. The parameters are obtained by matching theoretical vs empirical autocorrelation function. The theoretical autocorrelation function is computed according the methodology developed in Chadraa (2009).</p>
    </div>

    <div id="ref-usage">
    <div class="sourceCode"><pre class="sourceCode r"><code><span><span class="fu">gmm</span><span class="op">(</span><span class="va">yuima</span>, data <span class="op">=</span> <span class="cn">NULL</span>, <span class="va">start</span>,</span>
<span> method<span class="op">=</span><span class="st">"BFGS"</span>, fixed <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link">list</a></span><span class="op">(</span><span class="op">)</span>, <span class="va">lower</span>, <span class="va">upper</span>, lag.max <span class="op">=</span> <span class="cn">NULL</span>,</span>
<span> equally.spaced <span class="op">=</span> <span class="cn">FALSE</span>, aggregation<span class="op">=</span><span class="cn">TRUE</span>, Est.Incr <span class="op">=</span> <span class="st">"NoIncr"</span>, objFun <span class="op">=</span> <span class="st">"L2"</span><span class="op">)</span></span></code></pre></div>
    </div>

    <div id="arguments">
    <h2>Arguments</h2>
    <p></p>
<dl><dt id="arg-yuima">yuima<a class="anchor" aria-label="anchor" href="#arg-yuima"></a></dt>
<dd><p>a yuima object or an object of <code><a href="yuima.cogarch-class.html">yuima.cogarch-class</a></code>.</p></dd>

  <dt id="arg-data">data<a class="anchor" aria-label="anchor" href="#arg-data"></a></dt>
<dd><p>an object of class <code><a href="yuima.data-class.html">yuima.data-class</a></code> contains the observations available at uniformly spaced time. If <code>data=NULL</code>, the default, the function uses the data in an object of <code><a href="yuima-class.html">yuima-class</a></code>.</p></dd>

  <dt id="arg-start">start<a class="anchor" aria-label="anchor" href="#arg-start"></a></dt>
<dd><p>a <code>list</code> containing the starting values for the optimization routine.</p></dd>

  <dt id="arg-method">method<a class="anchor" aria-label="anchor" href="#arg-method"></a></dt>
<dd><p>a string indicating one of the methods available in <code><a href="https://rdrr.io/r/stats/optim.html" class="external-link">optim</a></code>.</p></dd>

  <dt id="arg-fixed">fixed<a class="anchor" aria-label="anchor" href="#arg-fixed"></a></dt>
<dd><p>a list of fixed parameters in optimization routine.</p></dd>

  <dt id="arg-lower">lower<a class="anchor" aria-label="anchor" href="#arg-lower"></a></dt>
<dd><p>a named list for specifying lower bounds of parameters.</p></dd>

  <dt id="arg-upper">upper<a class="anchor" aria-label="anchor" href="#arg-upper"></a></dt>
<dd><p>a named list for specifying upper bounds of parameters.</p></dd>

  <dt id="arg-lag-max">lag.max<a class="anchor" aria-label="anchor" href="#arg-lag-max"></a></dt>
<dd><p>maximum lag at which to calculate the theoretical and empirical acf. Default is <code>sqrt{N}</code> where <code>N</code> is the number of observation.</p></dd>

  <dt id="arg-equally-spaced">equally.spaced<a class="anchor" aria-label="anchor" href="#arg-equally-spaced"></a></dt>
<dd><p>Logical variable. If <code>equally.spaced = TRUE.</code>, the function use the returns of COGARCH(P,Q) evaluated at unitary length for the computation of the empirical autocorrelations. If <code>equally.spaced = FALSE</code>, the increments are evaluated on the interval with frequency specified in an object of class <code><a href="yuima.data-class.html">yuima.data-class</a></code> that contains the observed time series.</p></dd>

  <dt id="arg-aggregation">aggregation<a class="anchor" aria-label="anchor" href="#arg-aggregation"></a></dt>
<dd><p>If <code>aggregation=TRUE</code>, before the estimation of the levy parameters we aggregate the estimated increments</p></dd>

  <dt id="arg-est-incr">Est.Incr<a class="anchor" aria-label="anchor" href="#arg-est-incr"></a></dt>
<dd><p>a string variable, If <code>Est.Incr = "NoIncr"</code>, default value, <code>gmm</code> returns an object of class  <code><a href="cogarch.est.html">cogarch.est-class</a></code> that contains the COGARCH parameters.
  If <code>Est.Incr = "Incr"</code> or <code>Est.Incr = "IncrPar"</code> the output is an object of class <code>cogarch.est.incr-class</code>. In the first case the object contains the increments of underlying noise while in the second case also the estimated parameter of levy measure.</p></dd>

  <dt id="arg-objfun">objFun<a class="anchor" aria-label="anchor" href="#arg-objfun"></a></dt>
<dd><p>a string variable that identifies the objective function in the optimization step. <code>objFun = "L2"</code>, default value, the objective function is  a quadratic form where the weighting Matrix is the identity one. <code>objFun = "L2CUE"</code> the weighting matrix is estimated using Continuously Updating GMM (L2CUE).
  <code>objFun = "L1"</code>, the objective function is the mean absolute error. In the last case the standard error for estimators are not available.</p></dd>

</dl></div>
    <div id="details">
    <h2>Details</h2>
    <p>The routine is based on three steps: estimation of the COGARCH parameters, recovering the increments of the underlying Levy process and estimation of the levy measure parameters. The last two steps are available on request by the user.</p>
    </div>
    <div id="value">
    <h2>Value</h2>
    <p>The function returns a list with the same components of the object obtained when the function  <code><a href="https://rdrr.io/r/stats/optim.html" class="external-link">optim</a></code> is used.</p>
    </div>
    <div id="references">
    <h2>References</h2>
    <p>Chadraa, E. (2009) Statistical Modeling with COGARCH(P,Q) Processes. Phd Thesis</p>
    </div>
    <div id="author">
    <h2>Author</h2>
    <p>The YUIMA Project Team.</p>
    </div>

    <div id="ref-examples">
    <h2>Examples</h2>
    <div class="sourceCode"><pre class="sourceCode r"><code><span class="r-in"><span><span class="kw">if</span> <span class="op">(</span><span class="cn">FALSE</span><span class="op">)</span> <span class="op">{</span> <span class="co"># \dontrun{</span></span></span>
<span class="r-in"><span><span class="co"># Example COGARCH(1,1): the parameters are the same used in Haugh et al. 2005. In this case</span></span></span>
<span class="r-in"><span><span class="co"># we assume the underlying noise is a symmetric variance gamma.</span></span></span>
<span class="r-in"><span><span class="co"># As first step we define the COGARCH(1,1) in yuima:</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">mod1</span> <span class="op">&lt;-</span> <span class="fu"><a href="setCogarch.html">setCogarch</a></span><span class="op">(</span>p <span class="op">=</span> <span class="fl">1</span>, q <span class="op">=</span> <span class="fl">1</span>, work <span class="op">=</span> <span class="cn">FALSE</span>,</span></span>
<span class="r-in"><span>                   measure<span class="op">=</span><span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link">list</a></span><span class="op">(</span>df<span class="op">=</span><span class="st">"rbgamma(z,1,sqrt(2),1,sqrt(2))"</span><span class="op">)</span>,</span></span>
<span class="r-in"><span>                    measure.type <span class="op">=</span> <span class="st">"code"</span>, Cogarch.var <span class="op">=</span> <span class="st">"y"</span>,</span></span>
<span class="r-in"><span>                    V.var <span class="op">=</span> <span class="st">"v"</span>, Latent.var<span class="op">=</span><span class="st">"x"</span>,XinExpr<span class="op">=</span><span class="cn">TRUE</span><span class="op">)</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">param</span> <span class="op">&lt;-</span> <span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link">list</a></span><span class="op">(</span>a1 <span class="op">=</span> <span class="fl">0.038</span>,  b1 <span class="op">=</span>  <span class="fl">0.053</span>,</span></span>
<span class="r-in"><span>              a0 <span class="op">=</span> <span class="fl">0.04</span><span class="op">/</span><span class="fl">0.053</span>, x01 <span class="op">=</span> <span class="fl">20</span><span class="op">)</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="co"># We generate a trajectory</span></span></span>
<span class="r-in"><span><span class="va">samp</span> <span class="op">&lt;-</span> <span class="fu"><a href="setSampling.html">setSampling</a></span><span class="op">(</span>Terminal<span class="op">=</span><span class="fl">10000</span>, n<span class="op">=</span><span class="fl">100000</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="fu"><a href="https://rdrr.io/r/base/Random.html" class="external-link">set.seed</a></span><span class="op">(</span><span class="fl">210</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">sim1</span> <span class="op">&lt;-</span> <span class="fu"><a href="simulate.html">simulate</a></span><span class="op">(</span><span class="va">mod1</span>, sampling <span class="op">=</span> <span class="va">samp</span>, true.parameter <span class="op">=</span> <span class="va">param</span><span class="op">)</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="co"># We estimate the model</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">res1</span> <span class="op">&lt;-</span> <span class="fu">gmm</span><span class="op">(</span>yuima <span class="op">=</span> <span class="va">sim1</span>, start <span class="op">=</span> <span class="va">param</span><span class="op">)</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="fu"><a href="https://rdrr.io/r/base/summary.html" class="external-link">summary</a></span><span class="op">(</span><span class="va">res1</span><span class="op">)</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="op">}</span> <span class="co"># }</span></span></span>
</code></pre></div>
    </div>
  </div>
  <div class="col-md-3 hidden-xs hidden-sm" id="pkgdown-sidebar">
    <nav id="toc" data-toggle="toc" class="sticky-top"><h2 data-toc-skip>Contents</h2>
    </nav></div>
</div>


      <footer><div class="copyright">
  <p></p><p>Developed by YUIMA Project Team, Stefano M. Iacus.</p>
</div>

<div class="pkgdown">
  <p></p><p>Site built with <a href="https://pkgdown.r-lib.org/" class="external-link">pkgdown</a> 2.1.1.</p>
</div>

      </footer></div>






  </body></html>

