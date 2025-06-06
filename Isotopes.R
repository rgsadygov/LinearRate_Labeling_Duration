#
#
#
#   A method to compute isotope distribution of a 
#   peptide sequence given its amino acid sequence:
#   Elemental composition - maximum of 30 amino acids,
#   and maximum of 10 element types (C,H,N,O,S, X) for each element.
#
#    #1 Carbon, #2 Hydrogen, #3 Nitrogen, #4 Oxygen, #5 Sulfur
#


Letter2num <- function(Letter)
{
	utf8ToInt(Letter) - utf8ToInt("A") + 1L
}


ElementalComposition <- function(Sequence)
{


    szPeptide = as.character(Sequence);


    #print(c("seq = ", szPeptide));

    ElementalMatrix = matrix(0, nrow = 30, ncol=5);


    k = Letter2num("A");   

    ElementalMatrix[k, 1]  = 3;  ElementalMatrix[k, 2]  = 5;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    k = Letter2num("G");

    ElementalMatrix[k, 1]  = 2;  ElementalMatrix[k, 2]  = 3;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    k = Letter2num("S");

    ElementalMatrix[k, 1]  = 3;  ElementalMatrix[k, 2]  = 5;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 2;

    k = Letter2num("P");

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 7;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;


    k = Letter2num("V");

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 9;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;


    k = Letter2num("T");

    ElementalMatrix[k, 1]  = 4;  ElementalMatrix[k, 2]  = 7;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 2;


    k = Letter2num("C");

    ElementalMatrix[k, 1]  = 3;  ElementalMatrix[k, 2]  = 5;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;
    ElementalMatrix[k, 5]  = 1;


    k = Letter2num("L");

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 11;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;


    k = Letter2num("I");

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 11;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;


    k = Letter2num("N");

    ElementalMatrix[k, 1]  = 4;  ElementalMatrix[k, 2]  = 6;  ElementalMatrix[k, 3]  = 2;  ElementalMatrix[k, 4]  = 2;


    k = Letter2num("D");

    ElementalMatrix[k, 1]  = 4;  ElementalMatrix[k, 2]  = 5;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 3;


    k = Letter2num("Q");

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 8;  ElementalMatrix[k, 3]  = 2;  ElementalMatrix[k, 4]  = 2;


    k = Letter2num("K");

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 12;  ElementalMatrix[k, 3]  = 2;  ElementalMatrix[k, 4]  = 1;


    k = Letter2num("E");

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 7;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 3;


    k = Letter2num("M");

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 9;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    ElementalMatrix[k, 5]  = 1;


    k = Letter2num("H");

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 7;  ElementalMatrix[k, 3]  = 3;  ElementalMatrix[k, 4]  = 1;


    k = Letter2num("F");

    ElementalMatrix[k, 1]  = 9;  ElementalMatrix[k, 2]  = 9;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;


    k = Letter2num("R");

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 12;  ElementalMatrix[k, 3]  = 4;  ElementalMatrix[k, 4]  = 3;


    k = Letter2num("Y");

    ElementalMatrix[k, 1]  = 9;  ElementalMatrix[k, 2]  = 9;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 2;


    j = Letter2num("W");

    ElementalMatrix[j, 1]  = 11;  ElementalMatrix[j, 2]  = 10;  ElementalMatrix[j, 3]  = 2;  ElementalMatrix[j, 4]  = 1;



    nH = nO = nC = nN = nS = 0;

    nO = 1; nN = 0;

    nH = 3;
    

    for(i in 1:nchar(szPeptide))
    {
	 j = Letter2num(substring(szPeptide,i,i));

         nC = nC + as.numeric(ElementalMatrix[j, 1]);

	 nH = nH + ElementalMatrix[j, 2];
	 
         nN = nN + ElementalMatrix[j, 3];

	 nO = nO + ElementalMatrix[j, 4];

	 nS = nS + ElementalMatrix[j, 5];
    }

    Composition = list ();

    Composition$nC = nC;   Composition$nH = nH;
    Composition$nN = nN;   Composition$nO = nO;
    Composition$nS = nS;
    Composition$Sequence = szPeptide;
    Composition$Isotopes = Isotopes(nC, nH, nN, nO, nS);

    Composition
    
}


#
#
#   computes the isotope distributions given
#   the number of elemental atoms
#

Isotopes <- function(nC, nH, nN, nO, nS)
{

   #print(c(nC, nH, nN, nO, nS));
   
   probC = seq(1:32); probH = seq(1:32); probN = seq(1:32); probO = seq(1:32); probS = seq(1:32);

   for(i in 1:32)
     probC[i] = probH[i] = probN[i] = probO[i] = probS[i] = 0;

   pH = 0.00015574 / (0.99984426 + 0.00015574);

   probC[1] = 0.988922 / (0.988922 + 0.011078); probC[2] = 0.011078 / (0.988922 + 0.011078);


   probH[1] = 0.99984426 / (0.99984426 + 0.00015574); probH[2] = 0.00015574 / (0.99984426 + 0.00015574);


   probN[1] = 0.996337 / (0.996337 + 0.003663); probN[2] =  0.003663 / (0.996337 + 0.003663);


   probO[1] = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004); probO[2] = 0.000379 /(0.9976206 + 0.0003790 + 0.0020004);

   probO[3] = 0.0020004 /(0.9976206 + 0.0003790 + 0.0020004);


   probS[1] = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458); probS[2] = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

   probS[3] = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

   probS[5] = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);


   probTemp = probC;
   
   for(i in 1:(nC-1))
   {
	sum_all = sum(Re(fft(fft(probC)*fft(probTemp)/length(probC), inverse=TRUE)))

	probTemp = Re(fft(fft(probC)*fft(probTemp)/length(probC), inverse=TRUE)) / sum_all;    
   }

   probC = probTemp;


   probTemp = probH;
   
   for(i in 1:(nH-1))
   {
	sum_all = sum(Re(fft(fft(probH)*fft(probTemp)/length(probH), inverse=TRUE)))

	probTemp = Re(fft(fft(probH)*fft(probTemp)/length(probH), inverse=TRUE)) / sum_all;    
   }

   probH = probTemp;

   #print(c(probH[1:10], dbinom(0:9, nH, pH)));

   probTemp = probN;
   
   for(i in 1:(nN-1))
   {
	sum_all = sum(Re(fft(fft(probN)*fft(probTemp)/length(probN), inverse=TRUE)))

	probTemp = Re(fft(fft(probN)*fft(probTemp)/length(probN), inverse=TRUE)) / sum_all;    
   }

   probN = probTemp;


   probTemp = probO;
   
   for(i in 1:(nO-1))
   {
	sum_all = sum(Re(fft(fft(probO)*fft(probTemp)/length(probO), inverse=TRUE)))

	probTemp = Re(fft(fft(probO)*fft(probTemp)/length(probO), inverse=TRUE)) / sum_all;    
   }

   probO = probTemp;



   probTemp = probS;

   for(i in 1:(nS-1))
   {
	sum_all = sum(Re(fft(fft(probS)*fft(probTemp)/length(probS), inverse=TRUE)))

	probTemp = Re(fft(fft(probS)*fft(probTemp)/length(probS), inverse=TRUE)) / sum_all;    
   }

   sum_all = sum(Re(fft(fft(probC)*fft(probH)/length(probH), inverse=TRUE)))

   probTemp = Re(fft(fft(probC)*fft(probH)/length(probH), inverse=TRUE)) / sum_all;


   sum_all = sum(Re(fft(fft(probN)*fft(probTemp)/length(probN), inverse=TRUE)))

   probTemp = Re(fft(fft(probTemp)*fft(probN)/length(probN), inverse=TRUE)) / sum_all;    


   sum_all = sum(Re(fft(fft(probO)*fft(probTemp)/length(probO), inverse=TRUE)))

   probTemp = Re(fft(fft(probTemp)*fft(probO)/length(probO), inverse=TRUE)) / sum_all;    


   if(nS > 0)
   {
          sum_all = sum(Re(fft(fft(probS)*fft(probTemp)/length(probS), inverse=TRUE)))

	  probTemp = Re(fft(fft(probTemp)*fft(probS)/length(probS), inverse=TRUE)) / sum_all;    
   }


   probTemp
}


#
#
#     A function that computes isotope
#
#     distribution of a peptide sequence
#     to run it:
#
#     Isotope_Distribution("FVRRRAAAAAAAR")
#
#
#


Isotope_Distribution <- function(Sequence)
{
   
   Elements = ElementalComposition(Sequence);

   IsotopeDist = Isotopes(Elements$nC, Elements$nH, Elements$nN, Elements$nO,
   	       Elements$nS);


   # plot(c(0:(length(IsotopeDist)-1)), IsotopeDist, xlim=c(0, 15), xlab="Isotope Numbers", ylab="Relative Isotope Abundance",
   # 				      main=Sequence);

   # print(IsotopeDist);

}



#
#
#     A function that compares the factual and predicted
#
#     NEH and deuterium labeling using the new formula
#
#     to run it:
#
#     Check_Pxt("FVRRRAAAAAAAR", 0.0001, 30)
#
#     it should reproduce both: pxt  = 0.0001
#
#                         and:  NEH = 30;
#
#


Check_Pxt <- function(Sequence, pxt, NEH)
{
   
   pH = 0.00015574 / (0.99984426 + 0.00015574);

   pX = pH + pxt;

   probX = rep(0, 32);

   probX[1] = 1 - pX; probX[2] = pX;

   # The X element with pxt enrichment

   pTemp = probX;
   
   for(i in 1:(NEH-1))
   {
        sum_all = sum(Re(fft(fft(pTemp)*fft(probX)/length(probX), inverse=TRUE)))
	
   	pTemp = Re(fft(fft(pTemp)*fft(probX)/length(probX), inverse=TRUE)) / sum_all;  
   }


   probX = pTemp;

   print(probX[1:5]);

   xtemp = dbinom(0:(NEH-1), NEH, pX);

   print(xtemp[1:5]);

   #print(probX);

   Elements = ElementalComposition(Sequence);

   IsotopeDist = Isotopes(Elements$nC, (Elements$nH - NEH), Elements$nN, Elements$nO,
   	       Elements$nS);

   probTemp0 = Isotopes(Elements$nC, Elements$nH, Elements$nN, Elements$nO,
   	       Elements$nS);

   sum_all = sum(Re(fft(fft(probX)*fft(IsotopeDist)/length(probX), inverse=TRUE)))

   probTemp = Re(fft(fft(probX)*fft(IsotopeDist)/length(probX), inverse=TRUE)) / sum_all;  
   

   temp = (1 - pH)^2;

   #temp = temp * (probTemp[2]/probTemp[1] - probTemp0[2]/probTemp0[1]);

   #pxt_Computed = temp / (NEH + (1 - pH)*(probTemp[2]/probTemp[1] - probTemp0[2]/probTemp0[1]))

   temp = temp * (probTemp[2]*probTemp0[1] - probTemp0[2]*probTemp[1]);

   pxt_Computed = temp / (NEH*probTemp0[1]*probTemp[1] + (1 - pH)*(probTemp[2]*probTemp0[1] - probTemp0[2]*probTemp[1]))

   NEH_Computed = (1 - pH - pxt_Computed)*(probTemp[2]/probTemp[1] - probTemp0[2]/probTemp0[1])*(1 - pH)/pxt_Computed;

   print(c(pxt, pxt_Computed));

   print(c(NEH, NEH_Computed));

   print(c(probTemp[1:4], probTemp0[1:4]));

   A1_Ratio_t = probTemp[2]/probTemp[1];       # A_1(t) / A_0(t);

   A1_Ratio_0 = probTemp0[2]/probTemp0[1];     # A_1(0) / A_0(0); 

   A2_Ratio_t = probTemp[3] / probTemp[1];     #  = A_2(t)/A_0(t);

   A2_Ratio_0 = probTemp0[3] / probTemp0[1];   #  = A_2(0) / A_0(0);

   print(c(A2_Ratio_t, A2_Ratio_0, A1_Ratio_0, A1_Ratio_t));

   break;

   NFigure = 1000000;

   x_Figure = rep(110, NFigure);

   y_Figure = rep(110, NFigure);

   xmax = -1000; xmin=1000;

   imax = imin = 0;

   pxt_found = 0.0;

   x_Figure[1] = 0.0000001;

   for(i in 2:NFigure)
   {
       x_Figure[i] = x_Figure[i-1] + 0.0000001;


        NEH_Computed  = (1 - pH -x_Figure[i])*(A1_Ratio_t - A1_Ratio_0)*(1 - pH)/x_Figure[i];

	

       if( NEH_Computed  >  Elements$nH)
       {
          #next;
       }

      # print(i);

       y_Figure[i] = PxT_Formula(x_Figure[i], pH, A2_Ratio_t, A2_Ratio_0, A1_Ratio_0, A1_Ratio_t);

#       if((x_Figure[i] > pxt - 0.03*pxt) && (x_Figure[i] < pxt + 0.1*pxt))
       if((x_Figure[i] > pxt - 0.1*pxt) && (x_Figure[i] < pxt + 0.1*pxt))
       {

            print(c(x_Figure[i], y_Figure[i]));

	   if(x_Figure[i] > xmax)
	   {
		xmax = x_Figure[i];
		imax = i;
	   }

	   if(x_Figure[i] < xmin)
	   {
	        xmin = x_Figure[i];
		imin = i;
	    }


       }

       if(i >= 2 && (y_Figure[i] < 0.0 && y_Figure[i-1] > 0.0) ||
                 (y_Figure[i] > 0.0 && y_Figure[i-1] < 0.0) )
       {
            pxt_found = (x_Figure[i] + x_Figure[i-1])/2;

	    #x = uniroot(PxT_Formula, c(x_Figure[i-1], x_Figure[i+100]), tol = 10^(-17), pH = pH, A2_Ratio_t = A2_Ratio_t,
   	#		A2_Ratio_0 = A2_Ratio_0,  A1_Ratio_0 = A1_Ratio_0, A1_Ratio_t = A1_Ratio_t);

	    #print(x);
       }

       #print(c(x, x1));
   }

   print(c(length(x_Figure), length(y_Figure), max(y_Figure)));
   print(c("Found pxt = ", round(pxt_found, 7)));
   plot(x_Figure[imin:imax], y_Figure[imin:imax], xlab="Percent Body Water Enrichment",
   ylab = "Deviation from True Value");
#   lines(x_Figure[imin:imax], y_Figure[imin:imax], type="p");
   
  # print(c("Value for ", pxt, x1));

  # x = uniroot(PxT_Formula, c(0.0001, 0.1), tol = -10^(-9), pH = pH, NHtotal = Elements$nH, A2_Ratio_t = A2_Ratio_t,
   #			A2_Ratio_0 = A2_Ratio_0,  A1_Ratio_0 = A1_Ratio_0, A1_Ratio_t = A1_Ratio_t);

  # print(probTemp);

  # print(probTemp0);


   print(summary(lm(y_Figure[imin:imax] ~ x_Figure[imin:imax])));

   print(summary(lm(y_Figure[2:10000] ~ x_Figure[2:10000])));

#   print(x);
}

#
#
#   Function to compute the pxt from A2(t)/A0(t), A2(0)/A0(0), A1(t)/A0(t) only
#
#   no NEH will be used, the NEH value is actually passed via a formula
#
#
#

PxT_Formula <- function(px_t, pH, A2_Ratio_t, A2_Ratio_0, A1_Ratio_0, A1_Ratio_t)
{

    NEH = (1 - pH - px_t)*(A1_Ratio_t - A1_Ratio_0)*(1 - pH)/px_t;
#          (1 - pH - px_t)*(probTemp[2]/probTemp[1] - probTemp0[2]/probTemp0[1])*(1 - pH)/pxt;

    #print(c("NEH1 ", NEH[[1]]));


 

    f = -A2_Ratio_t + A2_Ratio_0 - A1_Ratio_0 * pH * NEH / (1 - pH);

    f = f - ((px_t + pH)/(1 - pH - px_t))^2 * NEH * (NEH + 1) / 2;

     f = f + (pH/(1 - pH))^2 * NEH * (NEH + 1)/2 ;

     f = f + A1_Ratio_t * (px_t + pH) * (A1_Ratio_t - A1_Ratio_0)*(1 - pH)/px_t ;

    f


}

Peptide_Pxt <- function(A2_Ratio_t, A2_Ratio_0, A1_Ratio_0, A1_Ratio_t, NEH)
{
   pH = 0.00015574 / (0.99984426 + 0.00015574);

   NEH_predicted = 0;

   pxt_predicted = 0.0;

   NFigure = 10000;


   x_Figure = rep(110, NFigure);

   y_Figure = rep(110, NFigure);

   pxt_found = 0.0; 

   beta_t = (A1_Ratio_t - A1_Ratio_0)/NEH

   c_t    = (1 - pH) * beta_t / (1 + (1 - pH)*beta_t);

   p_t    = (1 - pH) * c_t;

   print(c("Pxt from NEH: ", p_t));

   x_Figure[1] = 0.00001;

   y_Figure[1] = PxT_Formula(x_Figure[1], pH, A2_Ratio_t, A2_Ratio_0, A1_Ratio_0, A1_Ratio_t);

   for(i in 2:NFigure)
   {
       x_Figure[i] = x_Figure[i-1] + 0.00001;

       y_Figure[i] = PxT_Formula(x_Figure[i], pH, A2_Ratio_t, A2_Ratio_0, A1_Ratio_0, A1_Ratio_t);

       #print(c("i", i));


       NEH_temp = (1 - pH - x_Figure[i])*(A1_Ratio_t - A1_Ratio_0)*(1 - pH)/x_Figure[i];


       if((y_Figure[i] < 0.0 && y_Figure[i-1] > 0.0) ||
                 (y_Figure[i] > 0.0 && y_Figure[i-1] < 0.0) )
       {
            pxt_found = (x_Figure[i] + x_Figure[i-1])/2;

	    NEH_predicted = (1 - pH - x_Figure[i])*(A1_Ratio_t - A1_Ratio_0)*(1 - pH)/x_Figure[i];
           print(c("px and NEH: ", pxt_found, NEH));
       }

       if(NEH_temp > NEH-1 && NEH_temp < NEH+1)
       {
            #  print(c(x_Figure[i], y_Figure[i]));
       }
   }

   lines(x_Figure, y_Figure, type="p");

   print(c("pxt: NEH_predicted: ", pxt_found, NEH_predicted));

   print(lm(y_Figure ~ x_Figure));


   plot(x_Figure, y_Figure);
}