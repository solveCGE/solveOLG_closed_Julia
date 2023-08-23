# formatted reporting
function report(reporttext, reportcalc)
    cursorstart = 45
    
    countlinebreaks = 0
    for i in 1:length(reporttext)
        if reporttext[1:i] != "\n"
            break
        end
        countlinebreaks += 1
    end
    
    cursorstart = max(length(reporttext) - countlinebreaks, cursorstart) + 3 * (reportcalc >= 0) + 2 * (reportcalc < 0)
    
    charfill = join(fill(" ", cursorstart - length(reporttext) + countlinebreaks), "")
    
    println(reporttext, charfill, reportcalc)
end

# convert variable from cohort-view to period-view
function coh2per(inmat)
  maxage, numcoh = size(inmat)
  numper = numcoh-(maxage-1)
  
  if numper<=0
	error("coh2per: insufficient number of columns in input matrix")
  end
  
  outmat = zerosmat(maxage,numper)
  
  for a in 1:maxage
    outmat[a,:] = inmat[a,(maxage:numcoh).-(a-1)]
  end
  
  return outmat
end

# convert variable from cohort-view to period-view and aggregate over age
function aggcoh2per(inmat)
  return sum(coh2per(inmat),dims=1)
end

# convert variable from period-view to cohort-view
function per2coh(inmat)
  maxage, numper = size(inmat)
  numcoh = numper+(maxage-1)
  
  calibvec = inmat[:,1]
  
  outmat = zerosmat(maxage,numcoh)
  
  for a in 1:maxage
    if a < maxage
	  outmat[a,1:(maxage-a)] = onesrow(maxage-a)*calibvec[a]
	end
    outmat[a,(maxage:numcoh).-(a-1)] = inmat[a,:]
    if a > 1
	  outmat[a,(numcoh-(a-2)):numcoh] = onesrow(a-1)*inmat[a,numper]
	end
  end
  
  return outmat
end

function seq(from,to,lengthout=to-from+1)
	return collect(range(from,to,length=floor(Int,lengthout)))
end

function zeroscol(dim1)
	return zeros(dim1,1)
end

function zerosrow(dim1)
	return zeros(1,dim1)
end

function zerosmat(dim1,dim2)
	return zeros(dim1,dim2)
end

function onescol(dim1)
	return ones(dim1,1)
end

function onesrow(dim1)
	return ones(1,dim1)
end

function onesmat(dim1,dim2)
	return ones(dim1,dim2)
end