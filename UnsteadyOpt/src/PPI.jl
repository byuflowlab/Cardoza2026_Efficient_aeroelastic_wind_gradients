

function get_escalator(table, start_year::Int, start_month::Int, end_year::Int, end_month::Int)
    start_row = findfirst(x -> x == start_year, table["year"])
    start_row == nothing ? @warn("Start year $start_year not found in table.") : nothing

    end_row = findfirst(x -> x == end_year, table["year"])
    end_row == nothing ? @warn("End year $end_year not found in table.") : nothing

    cost_start = start_month == 13 ? table["annual"][start_row] : table["table"][start_row, start_month]
    cost_end = end_month == 13 ? table["annual"][end_row] : table["table"][end_row, end_month]
    return cost_end / cost_start
end

function get_GDP_escalator(start_year::Int, end_year::Int)
    start_row = findfirst(x -> x == start_year, GDP["year"])
    start_row == nothing ? @warn("Start year $start_year not found in table.") : nothing

    end_row = findfirst(x -> x == end_year, GDP["year"])
    end_row == nothing ? @warn("End year $end_year not found in table.") : nothing

    cost_start = GDP["absolute"][start_row]
    cost_end = GDP["absolute"][end_row]
    return cost_end / cost_start
end


#### PPI

# IPPI_FND
# self.escData["IPPI_FND"] = Escalator("Foundations", ["BHVY   "], [100.00]) #Each key might have multiple tables, in this case, just one table.
function IPPI_FND(start_year::Int, start_month::Int, end_year::Int, end_month::Int; scale=0.01)
    return 100 * get_escalator(BHVY, start_year, start_month, end_year, end_month) * scale
end

#IPPI_LPM
# self.escData["IPPI_LPM"] = Escalator("Permits, engineering (Land Based)", ["GDP"], [100.00])
function IPPI_LPM(start_year::Int, start_month::Int, end_year::Int, end_month::Int; scale=0.01)
    return 100 * get_GDP_escalator(start_year, end_year) * scale
end

#IPPI_LEL
# self.escData["IPPI_LEL"] = Escalator("Land Based Elect",["3353119", "335313P", "3359291  ", "GDP"], [40.00, 15.00, 35.00, 10.00]) #4 tables, with different weights.
function IPPI_LEL(start_year::Int, start_month::Int, end_year::Int, end_month::Int; scale=0.01)
    return (40 * get_escalator(T3353119, start_year, start_month, end_year, end_month) +
            15 * get_escalator(T335313P, start_year, start_month, end_year, end_month) +
            35 * get_escalator(T3359291, start_year, start_month, end_year, end_month) +
            10 * get_GDP_escalator(start_year, end_year)) * scale
end

#IPPI_RDC
# self.escData["IPPI_RDC"] = Escalator("Road & Civil Work", ["BHWY   "], [100.00])
function IPPI_RDC(start_year::Int, start_month::Int, end_year::Int, end_month::Int; scale=0.01)
    return 100 * get_escalator(BHWY, start_year, start_month, end_year, end_month) * scale
end

#IPPI_LAI
# self.escData["IPPI_LAI"] = Escalator("Land Based Assembly & installation", ["BHVY   "], [100.00])
function IPPI_LAI(start_year::Int, start_month::Int, end_year::Int, end_month::Int; scale=0.01)
    return 100 * get_escalator(BHVY, start_year, start_month, end_year, end_month) * scale
end

#IPPI_TPT
# self.escData["IPPI_TPT"] = Escalator("Transportation On/Offshore", ["4841212"], [100.00])
function IPPI_TPT(start_year::Int, start_month::Int, end_year::Int, end_month::Int; scale=0.01)
    return 100 * get_escalator(T4841212, start_year, start_month, end_year, end_month) * scale
end

#IPPI_LOM
# self.escData["IPPI_LOM"] = Escalator(["O&M Land Based                      ", ["GDP    "], [100.00]])
function IPPI_LOM(start_year::Int, start_month::Int, end_year::Int, end_month::Int; scale=0.01)
    return 100 * get_GDP_escalator(start_year, end_year) * scale
end


#IPPI_LLR
# self.escData["IPPI_LLR"] = Escalator(["Land Based Levelized Replacement    ", ["GDP    "], [100.00]])
function IPPI_LLR(start_year::Int, start_month::Int, end_year::Int, end_month::Int; scale=0.01)
    return 100 * get_GDP_escalator(start_year, end_year) * scale
end


#IPPI_LSE
# self.escData["IPPI_LSE"] = Escalator(["Land Based & Offshore Lease Cost    ", ["GDP    "], [100.00]])
function IPPI_LSE(start_year::Int, start_month::Int, end_year::Int, end_month::Int; scale=0.01)
    return 100 * get_GDP_escalator(start_year, end_year) * scale
end