#man braucht E_ikz, die Erwartungswerte fuer Item i, Summe k und Antwort in Kategorie z
#Category Probabilities for an item and a value of the latent variable
T_functions <- function(lavFit, RawData, Xi){
        RawData <- na.omit(RawData)
        Ests <- parameterEstimates(lavFit)
        T_Mat <- matrix(ncol=5)
        TempMat <- matrix(ncol=4)

        for (i in 1:ncol(RawData)){
                b <- get_thresholds(lavFit, names(RawData)[i])
                a <- get_slope(lavFit, names(RawData)[i])
                #b <- Ests[Ests$lhs==names(RawData[i]) & Ests$op=="|", names(Ests)=="est"]
                #a <- Ests[Ests$rhs==names(RawData[i]) & Ests$op=="=~",names(Ests)=="est"]

                #linearer Exponent fuer Kategorie c: z_c
                z_c <- function(b_c, a_i, Xi) {a_i*(Xi - b_c)}
                z <- z_c(b, a, Xi)

                #T-Funktionen (Kategorienwahrscheinlichkeiten fuer Item i und einen Wert von Xi
                T1 <- 1 - pnorm(z[[1]])
                T2 <- pnorm(z[[1]]) - pnorm(z[[2]])
                T3 <- pnorm(z[[2]]) - pnorm(z[[3]])
                T4 <- pnorm(z[[3]]) - pnorm(z[[4]])
                T5 <- pnorm(z[[4]])
                T_Mat <- rbind(T_Mat, c(T1, T2, T3, T4, T5))
                TempMat <- rbind(TempMat, b)
        }
        T_Mat <- T_Mat[2:nrow(T_Mat),]
        T_Mat
}


#sumscore probabilities for a value of the latent variable
sumscore_probs <- function(lavFit, RawData, Xi){
        numCat <- 5

        #vector that holds the sumscore-probabilities
        S <- c(1, rep(0, numCat*ncol(RawData)))

        #matrix that holds category probabilities (called T-funtions)
        CatProbs <- T_functions(lavFit, RawData, Xi)

        for (i in (1:ncol(RawData))){
                #matrix that will hold the results from the intermediate probability computations for the sumscores
                S_mat <- matrix(0, nrow=length(S), ncol=length(S))
                S_copy <- S
                for (s in i:(length(S)-numCat)){
                        for (k in 1:numCat){
                                S_copy[s+k] <-  S[s]*CatProbs[i,k]
                                S_mat[s+k, s] <- S_copy[s+k]

                        }
                }

                S <- rowSums(S_mat)
                S[i] <- 0
        }
        S

}


all_sumscore_probs <- function(lavFit, RawData){
        fscores <- predict(lavFit)
        fscores <- as.data.frame(fscores)
        Xi <- fscores[[1]]
        Mat <- matrix(ncol=51, byrow=TRUE)
        for (n in Xi){
                SProb <- sumscore_probs(lavFit, RawData, n)
                SProb <- c(SProb)
                Mat <- rbind(Mat, SProb)
        }
        Mat <- Mat[2:nrow(Mat),]
        Mat
}


#sumscore probabilities without item i for a value of the latent variable
#takes the column number of the item as input
sumscore_probs_without_item <- function(RawData, lavFit, Xi, Item){
        Data_noI <- na.omit(RawData)
        Data_noI <- Data_noI[c(-Item)]
        SnoI <- sumscore_probs(lavFit, Data_noI, Xi)
        return(SnoI)
}


#compute expected cell frequencies for one Item and one Scoregroup over all values of Xi
#Item needs to be a numeric value
expected_cell_frequencies <- function(RawData, lavFit, Item, Scoregroup){
        require(lavaan)

        k <- Scoregroup
        numCat <- 5
        #distribution of Xi
        fscores <- predict(lavFit)
        dist_xi <- table(round(fscores, 2))/length(fscores)
        dist_xi <- as.numeric(dist_xi)

        #matrix that will hold expected frequencies for the categories of the item and all values of Xi
        Res <- matrix(0, nrow=1, ncol=numCat, byrow=TRUE)

        for (xi in 1:length(dist_xi)){

                #compute category probabilities for all items and one value of Xi
                CatProbs <- T_functions(lavFit, RawData, dist_xi[xi])

                #compute sumscore probabilities without item
                S_star <- sumscore_probs_without_item(RawData, lavFit, dist_xi[xi], Item)

                #compute sumscore probabilities with all items
                S <- sumscore_probs(lavFit, RawData, dist_xi[xi])

                #vector that will hold expected cell frequencies for the item over all values of Xi
                E_z <- c()

                #compute expected value for answer in category z of the item when total score is k 
                for (z in 1:numCat){
                        numerator <- CatProbs[Item, z] * S_star[k - z] * dist_xi[xi]
                        #denominator <- S[k] * dist_xi[xi]
                        E_z <- c(E_z, numerator)  
                }
                if (length(E_z) < numCat) {E_z <- c(E_z, 0)}
                Res <- rbind(Res, E_z)
                Summe <- sum_cells(Res)
                print(sprintf("Iteration %s: Summe ist %s", xi, Summe))
        }
        #sum up cell frequencies for each value of Xi
        Res <- colSums(Res)
        Res

}


#compute expected cell frequencies for an item over all Scoregroups
all_expected_cell_frequencies  <- function(RawData, lavFit, Item){
        numCat <- 5
        Scoregroups <- 1:(numCat*ncol(RawData))
        Res <- matrix(0, nrow=length(Scoregroups), ncol=numCat, byrow=TRUE)

        #the Sumscore has to be at least equal to the number of items so we can start iterating there
        Start <- ncol(RawData)

        for (k in Start:length(Scoregroups)){
                exp_freqs <- expected_cell_frequencies(RawData, lavFit, Item, k)
                Res[k,] <- exp_freqs
        }

        rownames(Res) <- Scoregroups
        Res

}

#is this function still needed?
all_sumscore_probs_without_item <- function(lavFit, RawData, Item){
        fscores <- predict(lavFit)
        fscores <- as.data.frame(fscores)
        Xi <- fscores[[1]]
        Mat <- matrix(ncol=51, byrow=TRUE)

        for (n in Xi){
                SProb <- sumscore_probs_without_item(RawData, lavFit, n, Item)
                SProb <- c(n, SProb)
                Mat <- rbind(Mat, SProb)
        }
        Mat <- Mat[2:nrow(Mat),]
        return(Mat)
}


#observed freqencies for all items
all_observed_cell_frequencies <- function(RawData){
        require(reshape)

        Sumscore <- apply(na.omit(RawData), 1, sum)
        dat <- cbind(na.omit(RawData), Sumscore)
        dat.m <- melt(dat, id="Sumscore")
        observed <- ddply(dat.m, .variables="value", .fun=sum_cells)
        observed <- observed[order(observed$Sumscore, observed$variable),]

        observed
}

#compute observed freqs for one item
#takes raw data and item number as numeric value
observed_cell_frequencies <- function(RawData, Item){
        RawData <- na.omit(RawData)
        Item_name <- names(RawData[c(Item)])

        Formula <- sprintf("~ Sumscore + %s", Item_name)

        Sumscore <- apply(RawData, 1, sum)
        RawData <- cbind(RawData, Sumscore)

        Res <- xtabs(data=RawData, Formula)
        Res

}


#fill the observed cell frequencies into a 50 by 5 matrix
#takes Item number as numeric value
observed_freqs_matrix <- function(RawData, Item){
        numCat <- 5
        Obs <- observed_cell_frequencies(RawData, Item)
        Res <- matrix(0, nrow=numCat*ncol(RawData), ncol=numCat, byrow=TRUE)

        for (i in 1:nrow(Res)){
                for (s in 1:nrow(Obs)){
                        if (as.character(i) == rownames(Obs)[s]){
                                Res[i,] <- Obs[s, ]

                        }

                }

        }

        Res[is.na(Res)] <- 0
        Res

}


#collapsing the upper and lower rows
collapse <- function(Data){
        Head <- head(Data,5)
        Tail <- tail(Data,5)
        Head_new <- colSums(Head)
        Tail_new <- colSums(Tail)
        
        Res <- rbind(Head_new, Data[6:(nrow(Data)-5),], Tail_new)
        Res

}

#sum up all cells of a matrix
sum_cells <- function(Mat){
        Mat[is.na(Mat)] <- 0
        Summe <- 0
        for (i in 1:nrow(Mat)){
                for (j in 1:ncol(Mat)){
                        Summe <- Summe + Mat[i,j]

                }

        }
        Summe
}

#function to collapse adjacent rows when their sum is less than one
rec <- function(Mat){
        res <- c()
        ro <- 1

        while (ro < (nrow(Mat)-1)){

                if (sum(Mat[ro,]) < 1){
                        combineRows <- Mat[ro,] + Mat[ro+1,]
                        ro <- ro + 2
                        print("ro ist gelicht")
                        print(ro)

                        Mat <- Mat[ro:nrow(Mat), ]
                        Mat <- rbind(combineRows, Mat)
                        print("Mat ist:")
                        print(Mat)
                        res <- c(res, rec(Mat))
                        print("hat geklappt!")
                        print(res)
                }


                else {
                        print("Jo!")
                        res <- c(res, sum(Mat[ro,]))
                        ro <- ro + 1 
                        Mat <- Mat[ro:nrow(Mat),]
                        res <- c(res, rec(Mat))
                }
        }

        if (ro == ((nrow(Mat) -1)) & sum(Mat[ro+1,] < 1)){
                print("Wir sind da angelangt")
                combineRows <- Mat[ro,] + Mat[ro+1,]
                res <- c(res, sum(combineRows))
                return(res)


        }

        else if (ro == ((nrow(Mat) -1)) & sum(Mat[ro,] < 1)){
                combineRows <- Mat[ro,] + Mat[ro+1,]
                res <- c(res, sum(combineRows))
                return(res)


        }
        res


}

A <- matrix(c(0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4), ncol=3, byrow=TRUE)


sx2 <- function(RawData, Item, ExpMat, showIndMat=FALSE){
        numCat <- 5
        #read in raw data
        RawData <- na.omit(RawData)
        #compute matrix of observed cell frequencies
        O <- observed_freqs_matrix(RawData, Item)
        #assign zero to all NA-values
        O[is.na(O)] <- 0
        #this is the Matrix of expected cell frequencies
        E <- ExpMat
        E[is.na(E)] <- 0
        #get rid of first ten rows (because minimum sumscore has to be ten with 10 items and item categories form 1 to 5
        O <- O[10:nrow(O),]
        E <- E[10:nrow(E),]
        #collapse the first and last five rows, because expected scores ususally are two low
        O <- collapse(O)
        E <- collapse(E)
        E <- E*sum_cells(O)
        #naming after Kang & chen (2011): number of subjects with score k
        N_k <- apply(O, 1, sum)

        #computing a matrix of the same size E (and O)
        #but with a 1 where cell in E is below 1
        #and 0 where cell in E is above 1
        IndMat <- matrix(0, nrow=nrow(E), ncol=ncol(E), byrow=TRUE)
        for (ro in 1:nrow(IndMat)){
                for (co in 1:ncol(IndMat)){
                        if (E[ro,co] < 1) {IndMat[ro,co] <- 1}
                }
        }

        #now go through all rows from E
        #if sum of all cells below 1 is bigger than 1, sum them and start computing chi-squared statistic (obs - ex)**2/exp
        #else sum one more cell of the same row until sum is at least one and start computing chi-squared stat
        #IndMat -> IndVec
        SX2 <- c()
        Diag <- matrix(ncol=2, byrow=TRUE)
        for (ro in 1:nrow(IndMat)){
                IndVec <- c() #vector with cell indices of adjacent cells to be collapsed
                remaining <- c() #remaining cell indices that might also have to be collapsed to ensures minimum freq of 1

                #compute the two vectors
                for (co in 1:ncol(IndMat)){
                        if (IndMat[ro,co]==1){IndVec <- c(IndVec, co)}
                        else {remaining <- c(remaining, co)}

                }

                #compute SX2-statistic
                if (sum(E[ro,IndVec]) >= 1){
                        break
                }

                else {for (i in remaining){
                        IndVec <- c(IndVec, i)
                        IndMat[ro,i] <- 1
                        if (sum(E[ro,IndVec]) >= 1) {
                                break
                        }

                }

                }
                SX2 <- c(SX2 , (N_k[ro] * (sum(O[ro,IndVec]) - sum(E[ro,IndVec]))**2 / sum(E[ro,IndVec])))
                Diag <- rbind(Diag, c(sum(O[ro,IndVec]), sum(E[ro,IndVec])))
        }
        DF <- (numCat*ncol(RawData) - 2*numCat + 1) * numCat - 2 - sum_cells(IndMat)
        res <- c(SX2 = sum(SX2), DF = DF, numCollapsed = sum_cells(IndMat))
        print(SX2)
        print(Diag)
        if (showIndMat==TRUE) {print(IndMat)}
        res

}

#extract thresholds for an item from lavaan object
#takes fitted object and item name as string
get_thresholds <- function(lavFit, Item){
    Ests <- parameterEstimates(lavFit)
    thresh <- Ests[Ests$lhs==Item & Ests$op=="|", names(Ests)=="est"]

    thresh
}

#extract item slope parameter
#takes fitted object and item name as string
get_slope <- function(lavFit, Item){
    Ests <- parameterEstimates(lavFit)
    slope <- Ests[Ests$rhs==Item & Ests$op=="=~",names(Ests)=="est"]
    slope
}

#compute expected category frequencies for an item from a lavaan-object 
#takes the object and item name as string as inputs
expected_cat_freqs <- function(lavFit, Item){
    require(ggplot2)
    require(lavaan)

    #estimated values for the latent variable
    fscores <- predict(lavFit)
    fscores <- as.data.frame(fscores)
    Xi <- fscores[[1]]

    #estimated Parameters 
    DF1 <- parameterEstimates(lavFit)

    #get discrimination parameter from lavaan object
    #loading <- DF1[DF1$rhs==Item & DF1$op=="=~",names(DF1)=="est"]
    loading <- get_slope(lavFit, Item)

    #vector with thresholds
    thresh <- get_thresholds(lavFit, Item)

    #matrix that will hold category probabilities
    T_Mat <- matrix(ncol=5, byrow=TRUE)

    #linear exponent for category c
    z_c <- function(b_c, a, Theta) {a*(Theta - b_c)}

    #iterate over all values of Xi
    for (xi in Xi){

        #compute thresholds for the item and do this for each value of Xi
        z <- z_c(thresh, loading, xi)

        #T-functions: category probabilities for item i and a value of the latent variable
        T1 <- 1 - pnorm(z[[1]])
        T2 <- pnorm(z[[1]]) - pnorm(z[[2]])
        T3 <- pnorm(z[[2]]) - pnorm(z[[3]])
        T4 <- pnorm(z[[3]]) - pnorm(z[[4]])
        T5 <- pnorm(z[[4]])

        #collect all category probabilites in a matrix
        T_Mat <- rbind(T_Mat, c(T1, T2, T3, T4, T5))

    }
    T_Mat <- T_Mat[2:nrow(T_Mat),]

    MultinomData <- matrix(0, nrow=ncol(T_Mat), byrow=FALSE)

    for (ro in 1:nrow(T_Mat)){
        MultinomData <- cbind(MultinomData, rmultinom(1, 1, T_Mat[ro,]))

    }
    expectedCellCounts <- apply(MultinomData, 1, sum)
    expectedCellCounts

}

