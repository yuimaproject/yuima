<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Class for the Mathematical Description of Integral of a Stochastic Process — yuima.Integral-class • The YUIMA Project</title><!-- jquery --><script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js" integrity="sha512-v2CJ7UaYy4JwqLDIrZUI/4hqeoQieOmAZNXBeQyjo21dadnwR+8ZaIJVT8EE2iyI61OV8e6M8PP2/4hpQINQ/g==" crossorigin="anonymous" referrerpolicy="no-referrer"></script><!-- Bootstrap --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha256-bZLfwXAP04zRMK2BjiO8iu9pf4FbLqX6zitd+tIvLhE=" crossorigin="anonymous"><script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha256-nuL8/2cJ5NDSSwnKD8VqreErSWHtnEP9E7AySL+1ev4=" crossorigin="anonymous"></script><!-- bootstrap-toc --><link rel="stylesheet" href="../bootstrap-toc.css"><script src="../bootstrap-toc.js"></script><!-- Font Awesome icons --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/all.min.css" integrity="sha256-mmgLkCYLUQbXn0B1SRqzHar6dCnv9oZFPEC1g1cwlkk=" crossorigin="anonymous"><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/v4-shims.min.css" integrity="sha256-wZjR52fzng1pJHwx4aV2AO3yyTOXrcDW7jBpJtTwVxw=" crossorigin="anonymous"><!-- clipboard.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/clipboard.js/2.0.6/clipboard.min.js" integrity="sha256-inc5kl9MA1hkeYUt+EC3BhlIgyp/2jDIyBLS6k3UxPI=" crossorigin="anonymous"></script><!-- headroom.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/headroom.min.js" integrity="sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/jQuery.headroom.min.js" integrity="sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=" crossorigin="anonymous"></script><!-- pkgdown --><link href="../pkgdown.css" rel="stylesheet"><script src="../pkgdown.js"></script><meta property="og:title" content="Class for the Mathematical Description of Integral of a Stochastic Process — yuima.Integral-class"><meta property="og:description" content="The yuima.Integral class is a class of the  yuima package that extends the yuima-class it represents a integral of a stochastic process
 zt = int^{t}_0 h(theta, Xs, s) dXs"><!-- mathjax --><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" integrity="sha256-nvJJv9wWKEm88qvoQl9ekL2J+k/RWIsaSScxxlsrv8k=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/config/TeX-AMS-MML_HTMLorMML.js" integrity="sha256-84DKXVJXs0/F8OTMzX4UR909+jtl4G7SPypPavF+GfA=" crossorigin="anonymous"></script><!--[if lt IE 9]>
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
        <span class="version label label-default" data-toggle="tooltip" data-placement="bottom" title="">1.15.29</span>
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
    <h1>Class for the Mathematical Description of Integral of a Stochastic Process</h1>

    <div class="hidden name"><code>yuima.Integral-class.rd</code></div>
    </div>

    <div class="ref-description">
    <p>The <code>yuima.Integral</code> class is a class of the  <span class="pkg">yuima</span> package that extends the <code><a href="yuima-class.html">yuima-class</a></code> it represents a integral of a stochastic process</p>
<p><code> zt = int^{t}_0 h(theta, Xs, s) dXs</code></p>
    </div>


    <div id="slots">
    <h2>Slots</h2>
    <p>In the following we report the the additional slots of an object of class <code>yuima.Integral</code> with respect to the <code><a href="yuima-class.html">yuima-class</a></code>:</p><dl><dt><code>Integral</code>:</dt>
<dd><p>It is an object of class <code>Integral.sde</code> and it is composed by the following slots:</p><dl><dt><code>param.Integral</code>:</dt>
<dd><p>it is an object of class <code>param.Integral</code> and it is composed by the following slots:</p><dl><dt><code>allparam</code>:</dt>
<dd><p>labels of all parameters (model and  integral).</p></dd>

            <dt><code>common</code>:</dt>
<dd><p>common parameters.</p></dd>

            <dt><code>Integrandparam</code>:</dt>
<dd><p>labels of all parameters only in the integral.</p></dd>


</dl><p></p></dd>

        <dt><code>variable.Integral</code>:</dt>
<dd><p>it is an object of class <code>variable.Integral</code> and it is composed by the following slots:</p><dl><dt><code>var.dx</code>:</dt>
<dd><p>integral variable.</p></dd>

          <dt><code>lower.var</code>:</dt>
<dd><p>lower bound of support.</p></dd>

          <dt><code>upper.var</code>:</dt>
<dd><p>upper bound of support.</p></dd>

          <dt><code>out.var</code>:</dt>
<dd><p>labels of output.</p></dd>

          <dt><code>var.time</code>:</dt>
<dd><p>label of time.</p></dd>


</dl><p></p></dd>

        <dt><code>Integrand</code>:</dt>
<dd><p>it is an object of class <code>variable.Integral</code> and it is composed by the following slots:</p><dl><dt><code>IntegrandList</code>:</dt>
<dd><p>It is a <code>list</code> that contains the components of integrand <code>h(theta, Xs, s)</code>.</p></dd>

                <dt><code>dimIntegrand</code>:</dt>
<dd><p>a <code>numeric</code> object that is the dimensions of the output.</p></dd>


</dl><p></p></dd>


</dl><p></p></dd>


</dl></div>
    <div id="methods">
    <h2>Methods</h2>
    <p></p><dl><dt>simulate</dt>
<dd><p>simulation method. For more information see
	  <code><a href="simulate.html">simulate</a></code>.</p></dd>


</dl></div>

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

