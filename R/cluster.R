cluster <-
function(data)
{
  n <- 1
  for (i in 1 : (length(data) - 1))
    {
      if (data[i] != data[i + 1]) 
        {
          n[i + 1] <- n[i] + 1
        }
      else
        {
          n[i + 1] <- n[i]
        }
    }
  return(n)
}
