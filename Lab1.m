>> ((6.6-(3+3/14))*(5+5/6))/((21-1.25)/2.5)

ans =

    2.5000

>> help elfun
  Elementary math functions.
 
  Trigonometric.
    sin         - Sine.
    sind        - Sine of argument in degrees.
    sinh        - Hyperbolic sine.
    asin        - Inverse sine.
    asind       - Inverse sine, result in degrees.
    asinh       - Inverse hyperbolic sine.
    cos         - Cosine.
    cosd        - Cosine of argument in degrees.
    cosh        - Hyperbolic cosine.
    acos        - Inverse cosine.
    acosd       - Inverse cosine, result in degrees.
    acosh       - Inverse hyperbolic cosine.
    tan         - Tangent.
    tand        - Tangent of argument in degrees.
    tanh        - Hyperbolic tangent.
    atan        - Inverse tangent.
    atand       - Inverse tangent, result in degrees.
    atan2       - Four quadrant inverse tangent.
    atan2d      - Four quadrant inverse tangent, result in degrees.
    atanh       - Inverse hyperbolic tangent.
    sec         - Secant.
    secd        - Secant of argument in degrees.
    sech        - Hyperbolic secant.
    asec        - Inverse secant.
    asecd       - Inverse secant, result in degrees.
    asech       - Inverse hyperbolic secant.
    csc         - Cosecant.
    cscd        - Cosecant of argument in degrees.
    csch        - Hyperbolic cosecant.
    acsc        - Inverse cosecant.
    acscd       - Inverse cosecant, result in degrees.
    acsch       - Inverse hyperbolic cosecant.
    cot         - Cotangent.
    cotd        - Cotangent of argument in degrees.
    coth        - Hyperbolic cotangent.
    acot        - Inverse cotangent.
    acotd       - Inverse cotangent, result in degrees.
    acoth       - Inverse hyperbolic cotangent.
    hypot       - Square root of sum of squares.
    deg2rad     - Convert angles from degrees to radians.
    rad2deg     - Convert angles from radians to degrees.
 
  Exponential.
    exp         - Exponential.
    expm1       - Compute exp(x)-1 accurately.
    log         - Natural logarithm.
    log1p       - Compute log(1+x) accurately.
    log10       - Common (base 10) logarithm.
    log2        - Base 2 logarithm and dissect floating point number.
    pow2        - Base 2 power and scale floating point number.
    realpow     - Power that will error out on complex result.
    reallog     - Natural logarithm of real number.
    realsqrt    - Square root of number greater than or equal to zero.
    sqrt        - Square root.
    nthroot     - Real n-th root of real numbers.
    nextpow2    - Next higher power of 2.
 
  Complex.
    abs         - Absolute value.
    angle       - Phase angle.
    complex     - Construct complex data from real and imaginary parts.
    conj        - Complex conjugate.
    imag        - Complex imaginary part.
    real        - Complex real part.
    unwrap      - Unwrap phase angle.
    isreal      - True for real array.
    cplxpair    - Sort numbers into complex conjugate pairs.
 
  Rounding and remainder.
    fix         - Round towards zero.
    floor       - Round towards minus infinity.
    ceil        - Round towards plus infinity.
    round       - Round towards nearest integer.
    mod         - Modulus (signed remainder after division).
    rem         - Remainder after division.
    sign        - Signum.

>> help datafun
  Data analysis and Fourier transforms.
 
  Basic operations.
    max         - Largest component.
    min         - Smallest component.
    bounds      - Smallest and largest components.
    mean        - Average or mean value.
    median      - Median value.
    mode        - Mode, or most frequent value in a sample.
    std         - Standard deviation.
    var         - Variance.
    sort        - Sort in ascending order. 
    sortrows    - Sort rows in ascending order.
    issorted    - TRUE for sorted vector and matrices.
    sum         - Sum of elements.
    prod        - Product of elements.
    histogram   - Histogram.
    histcounts  - Histogram bin counts.
    trapz       - Trapezoidal numerical integration.
    movsum      - Moving sum of elements.
    movvar      - Moving variance.
    movstd      - Moving standard deviation.
    movmedian   - Moving median.
    movmean     - Moving mean.
    movmin      - Moving minimum.
    movmax      - Moving maximum.
    cumsum      - Cumulative sum of elements.
    cumprod     - Cumulative product of elements.
    cummin      - Cumulative smallest component.
    cummax      - Cumulative largest component.
    cumtrapz    - Cumulative trapezoidal numerical integration.
 
  Finite differences.
    diff        - Difference and approximate derivative.
    gradient    - Approximate gradient.
    del2        - Discrete Laplacian.
 
  Correlation.
    corrcoef    - Correlation coefficients.
    cov         - Covariance matrix.
    subspace    - Angle between subspaces.
 
  Filtering and convolution.
    filter      - One-dimensional digital filter.
    filter2     - Two-dimensional digital filter.
    conv        - Convolution and polynomial multiplication.
    conv2       - Two-dimensional convolution.
    convn       - N-dimensional convolution.
    deconv      - Deconvolution and polynomial division.
    detrend     - Linear trend removal.
 
  Fourier transforms.
    fft         - Discrete Fourier transform.
    fft2        - Two-dimensional discrete Fourier transform.
    fftn        - N-dimensional discrete Fourier Transform.
    ifft        - Inverse discrete Fourier transform.
    ifft2       - Two-dimensional inverse discrete Fourier transform.
    ifftn       - N-dimensional inverse discrete Fourier Transform.
    fftshift    - Shift zero-frequency component to center of spectrum.
    ifftshift   - Inverse FFTSHIFT.
    fftw        - Interface to FFTW library run-time algorithm tuning control.
 
  Missing data.
    ismissing          - Find missing data.
    standardizeMissing - Convert to standard missing data.
    rmmissing          - Remove standard missing data.
    fillmissing        - Fill standard missing data.
 
  Data preprocessing.
    isoutlier    - Find outliers in data.
    filloutliers - Replace outliers in data.

>> help graph2d
  Two dimensional graphs.
  
  Elementary X-Y graphs.
    plot      - Linear plot.
    loglog    - Log-log scale plot.
    semilogx  - Semi-log scale plot.
    semilogy  - Semi-log scale plot.
    polar     - Polar coordinate plot.
    plotyy    - Graphs with y tick labels on the left and right.
 
  Axis control.
    axis       - Control axis scaling and appearance.
    zoom       - Zoom in and out on a 2-D plot.
    grid       - Grid lines.
    box        - Axis box.
    rbbox      - Rubberband box.
    hold       - Hold current graph.
    axes       - Create axes in arbitrary positions.
    subplot    - Create axes in tiled positions.
 
  Graph annotation.
    plotedit  - Tools for editing and annotating plots.
    title     - Graph title.
    xlabel    - X-axis label.
    ylabel    - Y-axis label. 
    texlabel  - Produces the TeX format from a character string.
    text      - Text annotation.
    gtext     - Place text with mouse.
 
  Hardcopy and printing.
    print      - Print graph or Simulink system; or save graph to MATLAB file.
    printopt   - Printer defaults.
    orient     - Set paper orientation. 
 
  See also graph3d, specgraph.

>> help format
 format Set output format.
    format with no inputs sets the output format to the default appropriate
    for the class of the variable. For float variables, the default is
    format SHORT.
 
    format does not affect how MATLAB computations are done. Computations
    on float variables, namely single or double, are done in appropriate
    floating point precision, no matter how those variables are displayed. 
    Computations on integer variables are done natively in integer. Integer
    variables are always displayed to the appropriate number of digits for
    the class, for example, 3 digits to display the INT8 range -128:127.
    format SHORT and LONG do not affect the display of integer variables.
 
    format may be used to switch between different output display formats
    of all float variables as follows:
      format SHORT     Scaled fixed point format with 5 digits.
      format LONG      Scaled fixed point format with 15 digits for double
                       and 7 digits for single.
      format SHORTE    Floating point format with 5 digits.
      format LONGE     Floating point format with 15 digits for double and
                       7 digits for single.
      format SHORTG    Best of fixed or floating point format with 5 
                       digits.
      format LONGG     Best of fixed or floating point format with 15 
                       digits for double and 7 digits for single.
      format SHORTENG  Engineering format that has at least 5 digits
                       and a power that is a multiple of three
      format LONGENG   Engineering format that has exactly 16 significant
                       digits and a power that is a multiple of three.
 
    format may be used to switch between different output display formats
    of all numeric variables as follows:
      format HEX     Hexadecimal format.
      format +       The symbols +, - and blank are printed 
                     for positive, negative and zero elements.
                     Imaginary parts are ignored.
      format BANK    Fixed format for dollars and cents.
      format RAT     Approximation by ratio of small integers.  Numbers
                     with a large numerator or large denominator are
                     replaced by *.
 
    format may be used to affect the spacing in the display of all
    variables as follows:
      format COMPACT Suppresses extra line-feeds.
      format LOOSE   Puts the extra line-feeds back in.
 
    Example:
       format short, pi, single(pi)
    displays both double and single pi with 5 digits as 3.1416 while
       format long, pi, single(pi)
    displays pi as 3.141592653589793 and single(pi) as 3.1415927.
 
       format, intmax('uint64'), realmax
    shows these values as 18446744073709551615 and 1.7977e+308 while
       format hex, intmax('uint64'), realmax
    shows them as ffffffffffffffff and 7fefffffffffffff respectively.
    The HEX display corresponds to the internal representation of the value
    and is not the same as the hexadecimal notation in the C programming
    language.
 
    See also disp, display, isnumeric, isfloat, isinteger.

    Reference page for format

>> help disp
 disp Display array.
    disp(X) displays array X without printing the array name or 
    additional description information such as the size and class name.
    In all other ways it's the same as leaving the semicolon off an
    expression except that nothing is shown for empty arrays.
 
    If X is a string or character array, the text is displayed.
 
    See also fprintf, sprintf, int2str, num2str, rats, format, details.

    Reference page for disp
    Other functions named disp

>> disp MATLAB 
MATLAB

>> help sprintf
 sprintf Write formatted data to string or character vector
    STR = sprintf(FORMAT, A, ...) applies FORMAT to all elements of
    array A and any additional array arguments in column order, and returns
    the results as STR. FORMAT can be a character vector or a string
    scalar. The data type of STR is the same as the data type of FORMAT.
 
    [STR, ERRMSG] = sprintf(FORMAT, A, ...) returns an error message when
    the operation is unsuccessful.  Otherwise, ERRMSG is empty.
 
    sprintf is the same as FPRINTF except that it returns the data in a 
    MATLAB variable rather than writing to a file.
 
    FORMAT describes the format of the output fields, and can include 
    combinations of the following:
 
       * Conversion specifications, which include a % character, a
         conversion character (such as d, i, o, u, x, f, e, g, c, or s),
         and optional flags, width, and precision fields.  For more
         details, type "doc sprintf" at the command prompt.
 
       * Literal text to print.
 
       * Escape characters, including:
             \b     Backspace            ''   Single quotation mark
             \f     Form feed            %%   Percent character
             \n     New line             \\   Backslash
             \r     Carriage return      \xN  Hexadecimal number N
             \t     Horizontal tab       \N   Octal number N%
         where \n is a line termination character on all platforms.
 
    Notes:
 
    If you apply an integer or text conversion to a numeric value that
    contains a decimal, MATLAB overrides the specified conversion, and
    uses %e to express the value in exponential notation.
 
    Numeric conversions print only the real component of complex numbers.
 
    Examples
       sprintf('%0.5g',(1+sqrt(5))/2)       % 1.618
       sprintf('%0.5g',1/eps)               % 4.5036e+15       
       sprintf('%15.5f',1/eps)              % 4503599627370496.00000
       sprintf('%d',round(pi))              % 3
       sprintf('%s','hello')                % hello
       sprintf('The array is %dx%d.',2,3)   % The array is 2x3.
 
    See also fprintf, sscanf, num2str, int2str, char, string, compose.

    Reference page for sprintf
    Other functions named sprintf

>> r=1

r =

     1

>> sprintf('Радіус r=%4.2f м\n', r)

ans =

    'Радіус r=1.00 м
     '

>> sprintf('Довжина окружності L=%4.3f м\n', 2*pi*r)

ans =

    'Довжина окружності L=6.283 м
     '

>> fprintf('Радіус r=%4.2f м\n', r)
Радіус r=1.00 м
>> fprintf('Довжина окружності L=%4.3f м\n', 2*pi*r)
Довжина окружності L=6.283 м
>> r=[1 2 3];
>> fprintf('Радіус r=%4.2f м\n', r)
Радіус r=1.00 м
Радіус r=2.00 м
Радіус r=3.00 м
>> fprintf('Довжина окружності L=%4.3f м\n', 2*pi*r)
Довжина окружності L=6.283 м
Довжина окружності L=12.566 м
Довжина окружності L=18.850 м
>> a=4.13*(10^(-1))

a =

    0.4130

>> b=1/261

b =

    0.0038

>> a=1.5

a =

    1.5000

>> b=0.8

b =

    0.8000

>> A=61

A =

    61

Q= sqrt((a*(sqrt(b)))/((tan (A))^(1/3)))

Q =

    0.9296