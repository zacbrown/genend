#!/usr/bin/env ruby

require 'rubygems'
require 'gruff'
#require 'gsl'

SPACINGS = 50

class Record
  attr_accessor :score, :match

  def initialize(score, match)
    @score = score
    @match = match
  end
end

class Array
  def count_true
    ret_val = 0
    self.each { |val|
      ret_val += 1 if val == true
    }
    return ret_val
  end
end

def process_groups(score_array, match_array)
  small = score_array[0]
  len = score_array.size
  ret_array = [[],[],[],[],[],[]]
  match_index = 0

  match_array.each { |row|
    val = row.index('true')
    if val == nil
      ret_array.each { |array| array << small }
    else
      ret_array[val] << score_array[match_index]
      match_index += 1

      for i in (0..5)
        if i == val then next end
        ret_array[i] << small
      end
    end
  }

  return ret_array
end

# def poly_fit(data, degree)
#   y = Vector(data)
#   last_val = data.length - 1
#   x = Vector((0..last_val).to_a)

#   coef, err, chisq, status = MultiFit.polyfit(x, y, degree)
#   x2 = Vector.linspace(1, last_val, SPACINGS)

#   return coef.eval(x2)
# end

file = ARGV[0]
split_num = ARGV[1]
piece_size = file.split('-')[1]
date = file.split('-')[2].chomp('.dat')
header = file.split('-')[0]
NUM_KMERS = 31050

score_vals = {'3'=>[],'4'=>[],'5'=>[],'6'=>[],'7'=>[],'8'=>[]}

File.open(file, 'r').each { |line|
  tokens = line.split(' ')
  score_vals[tokens[0]] << Record.new(tokens[7].to_f, tokens[1..6])
}

graph_hash = {}
low = 0.0; high = 0.0;
short_len = 1000000000

score_vals.each { |key,value|
  next if value == []
#  value.sort! { |a,b| a.score <=> b.score }
  puts "processing: #{key}-mers"

  score_arr = Array.new
  match_arr = Array.new
  value.each { |record|
    score_arr << record.score
    match_arr << record.match
  }

  short_len = score_arr.size if score_arr.size < short_len
#   low = ("%5.2f" % score_arr[0]).to_f if key == '3'
#   high = ("%5.2f" % score_arr[score_arr.size - 2]).to_f if key == '8'
  graph_hash[key] = process_groups(score_arr, match_arr)
}

# truncate arrays that are longer than the shortest in the group
# graph_hash.each { |key,value|
#   if value.size > short_len
#     diff = value.size - short_len
#     value.slice! 0, diff
#   end
# }

colors = {:my_green =>'#00FF00',:my_pale_green => '#66CC00',
  :my_pale_yellow => '#CCFF00', :my_yellow => '#FFFF00',
  :my_orange => '#FFA500', :my_red => '#FF0000'}

# setup data
graph_hash.each { |key,value|
  graph = Gruff::Area.new('1200x800')
  graph.theme_keynote #2E37FE
  graph.replace_colors ["#ffffff", "#00ff00", "#3333FF",
                      "#ff00ff", "#ffff00", "#ff0000"]
#   graph.replace_colors [colors[:my_red],
#                         colors[:my_orange],
#                         colors[:my_yellow],
#                         colors[:my_pale_yellow],
#                         colors[:my_pale_green],
#                         colors[:my_green]]
  graph.title = "#{header} #{key}-mer 10000 bp Overlap Taxonomic Classification"
  graph.title_font_size = 20
  graph.marker_font_size = 16
  graph.x_axis_label = "relative piece location"
  graph.y_axis_label = "raw score"
  # graph.maximum_value = 100
  # graph.minimum_value = 0

  array_dat = value.reverse
  graph.data "phylum", array_dat[0]
  graph.data "class", array_dat[1]
  graph.data "order", array_dat[2]
  graph.data "family", array_dat[3]
  graph.data "genus", array_dat[4]
  graph.data "species", array_dat[5]

  graph.write("result-#{header}-#{piece_size}-#{key}-#{date}.png")

}


