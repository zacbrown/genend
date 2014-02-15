#!/usr/bin/env ruby

require 'rubygems'
gem PLATFORM == 'java' ? 'rmagick4j' : 'rmagick'
require 'RMagick'
#require 'gsl'

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

  def sum
    inject( nil ) { |sum,x| sum ? sum+x : x }
  end

  def mean
    sum / size
  end
end

def process_groups(score_array, match_array)
  small = score_array.mean
  len = score_array.size
  ret_array = [[],[]]
  arrindices = ["s","g","f","o","c","p"]
  match_index = 0

  match_array.each { |row|
    ind = 0
    row.each do |val|
      if val == "true"
        ret_array[0] << arrindices[ind]
        break
      end
      ind += 1
    end
    if ind == 6
      ret_array[0] << "n"
    end
    ret_array[1] << score_array[match_index]
    match_index += 1
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
spec_names = {
  "NC_009613"=>"Flavobacterium psychrophilum",
  "NC_004603"=>"Vibrio parahaemolyticus",
  "NC_003295"=>"Ralstonia solanacearum",
  "NC_007622"=>"Staphylococcus aureus",
  "NC_005823"=>"Leptospira interrogans"
}

File.open(file, 'r').each { |line|
  tokens = line.split(' ')
  score_vals[tokens[0]] << Record.new(tokens[7].to_f, tokens[1..6])
}

graph_hash = {}
short_len = 1000000000

score_vals.each { |key,value|
  next if value == []
  puts "processing: #{key}-mers"

  score_arr = Array.new
  match_arr = Array.new
  value.each { |record|
    score_arr << record.score
    match_arr << record.match
  }

  short_len = score_arr.size if score_arr.size < short_len
  graph_hash[key] = process_groups(score_arr, match_arr)
}

colors = {:my_green =>'#00FF00',:my_pale_green => '#66CC00',
  :my_pale_yellow => '#CCFF00', :my_yellow => '#FFFF00',
  :my_orange => '#FFA500', :my_red => '#FF0000'}

# ["#ffffff", "#00ff00", "#3333FF",
 #                     "#ff00ff", "#ffff00", "#ff0000"]

graph_hash.each do |key, value|
  canvas = Magick::Image.new(4000, 800,
                             Magick::HatchFill.new('white','lightcyan2'))
  gc = Magick::Draw.new
  gc.stroke_width(2)
  x0 = 75
  y0 = 750
  x = x0
  gc.line(x0, y0, x0+3900, y0)
  gc.line(x0, y0, x0, y0-700)

  match = value[0]
  score = value[1]
  ind = 0
  match.each do |row|
    case row
      when "s"
      gc.stroke('#ff0000')
      gc.fill('#ff0000')
      when "g"
      gc.stroke('#ffff00')
      gc.fill('#ffff00')
      when "f"
      gc.stroke('#ff00ff')
      gc.fill('#ff00ff')
      when "o"
      gc.stroke('#3333ff')
      gc.fill('#3333ff')
      when "c"
      gc.stroke('#00ff00')
      gc.fill('#00ff00')
      when "p"
      gc.stroke('#8B8B00')
      gc.fill('#8B8B00')
      else
      gc.stroke('#000000')
      gc.fill('#000000')
    end
    gc.fill_opacity(1)
    norm_score = ("%5.6f" % ((score[ind] /= 10000)+3)).to_f
    y = y0 - (norm_score * 200)
    gc.circle(x, y, x+2, y)
    x += 5
    ind += 1
  end

  #draw legend
  gc.fill_opacity(1)
  gc.fill('#ff0000')
  gc.circle(300, 100, 305, 100)
  gc.stroke('transparent').fill('black')
  gc.text(315, 100, "species match")

  gc.fill_opacity(1)
  gc.fill('#ffff00')
  gc.circle(300, 115, 305, 115)
  gc.stroke('transparent').fill('black')
  gc.text(315, 115, "genus match")

  gc.fill_opacity(1)
  gc.fill('#ff00ff')
  gc.circle(300, 130, 305, 130)
  gc.stroke('transparent').fill('black')
  gc.text(315, 130, "family match")

  gc.fill_opacity(1)
  gc.fill('#3333ff')
  gc.circle(300, 145, 305, 145)
  gc.stroke('transparent').fill('black')
  gc.text(315, 145, "order match")

  gc.fill_opacity(1)
  gc.fill('#00ff00')
  gc.circle(300, 160, 305, 160)
  gc.stroke('transparent').fill('black')
  gc.text(315, 160, "class match")

  gc.fill_opacity(1)
  gc.fill('#8B8B00')
  gc.circle(300, 175, 305, 175)
  gc.stroke('transparent').fill('black')
  gc.text(315, 175, "phylum match")

  gc.fill_opacity(1)
  gc.fill('#000000')
  gc.circle(300, 190, 305, 190)
  gc.stroke('transparent').fill('black')
  gc.text(315, 190, "unmatched")

  #draw title
  gc.stroke('transparent').fill('black')
  gc.pointsize(48)
  gc.text(750, 190, "#{spec_names[header]} 10000bp #{key}-mers #{date}")

  gc.draw(canvas)
  canvas.write("result-#{header}-10000-#{key}-#{date}.png")
end
